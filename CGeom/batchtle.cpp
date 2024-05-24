/*
	batchtle.cpp
	combined
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <omp.h>
#include <mutex>
#include "CGcyl.h"
#include "CGrpp.h"
#include "CGtrc.h"
#include "CGset.h"
#include "UVfile.h"
#include "Geometry.h"
#include "TRcoord.h"
#include "Scatter.h"
#include "Tally.h"
#include "Xor128.h"

int main()
{
  double wtime0 = omp_get_wtime();
  int prlev = 0;
  Geometry *geom = new Geometry();
  /*  CG description  */
  int ncgs = 3;
  CG **cg = new CG*[ncgs];
  cg[0] = new CGrpp();
  double bmin[] = {-17.55,-17.55,-0.1};
  double bmax[] = { 17.55, 17.55,10.56};
  cg[0]->RPP(bmin,bmax);
  cg[1] = new CGtrc();
  cg[1]->TRC(8.0,17.55,1.0,9.55);
  cg[2] = new CGcyl();
  cg[2]->CYL(8.0, 0.0,1.0);
  geom->setCG(ncgs,cg);

  /*  UV description  */
  int nxyz[] = {128,128,134};
  int nvoxel = nxyz[0]*nxyz[1]*nxyz[2];
  double dxyz[] = {0.2,0.2,0.2};
  UVfile *uvf = new UVfile();
  uvf->open("/home/jhchang/Dicom3/test3/0831150844478.uv",nxyz,dxyz);
  geom->setLattice(nxyz,dxyz,uvf->UV());

  TRcoord *trc = new TRcoord();
  double Zb = 20.0; // 20.0;
  double theta = 85.0; // 45.0;
  double phi = 0.0;
  double psi = 0.0;
  double target[] = {12.0,12.0,25.0};
  trc->set(target,Zb,theta,phi,psi);

  /*  set source coordinate  */
  double CGp[] = {0.0,0.0,-1.0};
  double CGw[] = {0.1,-0.1,1.0};
  CGw[2] = sqrt(1.0-CGw[0]*CGw[0]-CGw[1]*CGw[1]);
  double Rsrc = 8.0;

  Scatter *scat = new Scatter();

  Tally *tally = new Tally();
  tally->setLattice(nxyz,dxyz,uvf->UV());

  init_xor128();
  for(int nj=0; nj < 7; nj++) jump_xor128();

  double PI = 2*acos(0.0);
  int nbatch = 2000; // 2000;	// batch size
  int nhist = 200;	// no. of batchs
  /*  total tally  */
  float *tottal = new float[nvoxel];
  for(int ijk=0; ijk < nvoxel; ijk++) tottal[ijk] = 0.0;
  // std::mutex tally_mutex;
  int cpu = 0;
  double cpusum = 0;
  double cpu2sum = 0;
  #pragma omp firstgprivate(geom,trc,tally, CGp,CGw)
  {
  #pragma omp parallel for
  for(int nh=0; nh < nhist; nh++) {

    double UVp[3],UVw[3];

    clock_t start = clock();
    tally->clear();
    for(int iter=0; iter < nbatch; iter++) {
      /*  generate source position in disk */
      double xsrc,ysrc;
      double rsq = 99.0;
      while(rsq > 1.0) {
        xsrc = 2*ranf() - 1.0;
        ysrc = 2*ranf() - 1.0;
        rsq = xsrc*xsrc+ysrc*ysrc;
      }
      CGp[0] = Rsrc*xsrc;  CGp[1] = Rsrc*ysrc;
      /*  generate azimuthal angle  */
      double sintheta = sqrt(1.0-CGw[2]*CGw[2]);
      double ang = 2*PI*ranf();
      CGw[0] = sintheta*cos(ang);   CGw[1] = sintheta*sin(ang);

      trc->S2U(CGp,UVp);
      trc->S2Uvec(CGw,UVw);
      geom->setTR(trc);
      int old = -1;
      double eps = 1.0e-10;
      int ncross = 0;
      geom->setPrlev(0);
      int neps = 10;
      while(1) {
        Ray ray = geom->track(UVp,UVw);
        int isuv = 1;
        if(ray.here < 32) isuv = 0;

        if(ray.here < 0) break;
        double dist = ray.dist;
        double fflight = 1.0e+20;
        int vacuum = 0;
        if(ray.here == 0) vacuum = 1;
        else if(ray.here == 32) vacuum = 1;
        if(!vacuum) {
          /*  compute flight distance  -ln(xsi)/Sig_t  */
          double sigtot = 100.0;
          if(ray.here < 32) sigtot = 0.2;
          fflight = -log(ranf())/sigtot;
        }
 
        int iso = 0;
        if(fflight < dist) {
          dist = fflight;
          ray.next = ray.here;
          if(prlev) printf("# collide dist=%.2le here=%d(%c)\n",dist,ray.here,ray.here);
          iso = 1;  // isotropic
        }
        if(dist < eps) dist = eps;
        else {
          if(isuv && !vacuum) tally->TLE(UVp,UVw, ray.here, dist);
        }
        /*  change position */
        for(int j=0; j < 3; j++) UVp[j] += dist*UVw[j];
        /*  change direction */
        if(iso) {
          scat->isotropic(1.0,1.0,UVw);
        }
        if(ray.next < 0) {
          break;
        }
        old = ray.next;
        ncross++;
      }	// while 1

    }	// for iter
    // tally_mutex.lock();
    // #pragma omp critical
    for(int ijk=0; ijk < nvoxel; ijk++) tottal[ijk] += tally->tal[ijk];
    // tally_mutex.unlock();
    double cpu_time = (double) (clock() - start) / CLOCKS_PER_SEC;
    printf("%d: %.4lf\n",nh, cpu_time);
    cpu++;
    cpusum += cpu_time;
    cpu2sum += cpu_time*cpu_time;
  }	// for nh
  }
  printf("=== done ===\n");
  tally->print(tottal);
  printf(" nbatch=%d batchs=%d: ",nbatch,nhist);
  double cpuavg = cpusum/cpu;
  printf("  cpu avg=%.4lf sec",cpuavg);
  double cpuvar = cpu2sum/cpu - cpuavg*cpuavg;
  printf("  var=%.4lf  std=%.2lf %%\n",cpuvar, 100*sqrt(cpuvar)/cpuavg);
  /*  find total track length */
  float tracktot = 0;
  float trackmax = 0;
  for(int ijk=0; ijk < nvoxel; ijk++) {
    float t = tottal[ijk];
    if(t > trackmax) trackmax = t;
    tracktot += t;
  }
  printf(" total track length=%.3le max=%.3le\n",tracktot,trackmax);
  printf(" wall time=%.6lf\n", omp_get_wtime()-wtime0);


  delete [] tottal;
  delete tally;
  delete trc;
  delete scat;
  delete uvf;
  delete geom;
  for(int n=0; n < ncgs; n++) delete cg[n];
  delete [] cg;
};
