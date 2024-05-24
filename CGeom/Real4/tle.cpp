/*
	tle.cpp - float version
	combined
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "CGcyl.h"
#include "CGrpp.h"
#include "CGtrc.h"
#include "CGset.h"
#include "UVfile.h"
#include "Geometry.h"
#include "TRcoord.h"
#include "Tally.h"
#include "Xor128.h"

int main()
{
  int prlev = 0;
  Geometry *geom = new Geometry();
  /*  CG description  */
  int ncgs = 3;
  CG **cg = new CG*[ncgs];
  cg[0] = new CGrpp();
  float bmin[] = {-17.55,-17.55,-0.1};
  float bmax[] = { 17.55, 17.55,10.56};
  cg[0]->RPP(bmin,bmax);
  cg[1] = new CGtrc();
  cg[1]->TRC(8.0,17.55,1.0,9.55);
  cg[2] = new CGcyl();
  cg[2]->CYL(8.0, 0.0,1.0);
  geom->setCG(ncgs,cg);

  /*  UV description  */
  int nxyz[] = {128,128,134};
  float dxyz[] = {0.2,0.2,0.2};
  UVfile *uvf = new UVfile();
  uvf->open("/home/jhchang/Dicom3/test3/0831150844478.uv",nxyz,dxyz);
  geom->setLattice(nxyz,dxyz,uvf->UV());

  TRcoord *trc = new TRcoord();
  float Zb = 20.0; // 20.0;
  float theta = 85.0; // 45.0;
  float phi = 0.0;
  float psi = 0.0;
  float target[] = {12.0,12.0,25.0};
  trc->set(target,Zb,theta,phi,psi);

  /*  set source coordinate  */
  float CGp[] = {0.0,0.0,-1.0};
  float CGw[] = {0.1,-0.1,1.0};
  CGw[2] = sqrt(1.0-CGw[0]*CGw[0]-CGw[1]*CGw[1]);

  float UVp[3],UVw[3];
  CGp[0] = -20.0;   CGp[1] = CGp[0];
  Tally *tally = new Tally();
  tally->setLattice(nxyz,dxyz,uvf->UV());
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();

  for(int iter=0; iter < 1000; iter++) {
    if(prlev > 1) {
      printf("### CGp=(%.3lf,%.3lf,%.3lf) ",CGp[0],CGp[1],CGp[2]);
      printf("    CGw=(%.3lf,%.3lf,%.3lf) ",CGw[0],CGw[1],CGw[2]);
      printf("\n");
    }
    trc->S2U(CGp,UVp);
    trc->S2Uvec(CGw,UVw);
    geom->setTR(trc);
    int old = -1;
    float eps = 1.0e-5;	// for float precision
    int ncross = 0;
    geom->setPrlev(0);
    int neps = 10;
    while(1) {
      // printf("=== UVp=(%.3lf,%.3lf,%.3lf) ",UVp[0],UVp[1],UVp[2]);
      // printf("   UVw=(%.3lf,%.3lf,%.3lf)\n",UVw[0],UVw[1],UVw[2]);
      Ray ray = geom->track(UVp,UVw);
      int isuv = 1;
      if(ray.here < 32) isuv = 0;

      if(ray.here < 0) {
        if(prlev) printf("### escape\n");
        break;
      }
      float dist = ray.dist;
      float fflight = 1.0e+20;
      int vacuum = 0;
      if(ray.here == 0) vacuum = 1;
      else if(ray.here == 32) vacuum = 1;
      if(!vacuum) {
        /*  compute flight distance  -ln(xsi)/Sig_t  */
        float sigtot = 100.0;
        if(ray.here < 32) sigtot = 0.2;
        fflight = -log(ran->ranf())/sigtot;
      }
 
      int iso = 0;
      if(fflight < dist) {
        dist = fflight;
        ray.next = ray.here;
        if(prlev) printf("# collide dist=%.2le here=%d(%c)\n",dist,ray.here,ray.here);
        iso = 1;  // isotropic
      }
      if(prlev) {
        printf("%d %.3lf %.3lf %.3lf  %.2lf",ray.here,UVp[0],UVp[1],UVp[2], dist);
        printf("  %.3lf %.3lf %.3lf",dist*UVw[0],dist*UVw[1],dist*UVw[2]);
        printf("  %.3lf %.3lf %.3lf",UVp[0]+dist*UVw[0],UVp[1]+dist*UVw[1],UVp[2]+dist*UVw[2]);
        printf("\n");
      }
      if(dist < eps) dist = eps;
      else {
        if(isuv && !vacuum) tally->TLE(UVp,UVw, ray.here, dist);
      }
      /*  change position */
      for(int j=0; j < 3; j++) UVp[j] += dist*UVw[j];
      /*  change direction */
      if(iso) {
        ran->direcS2(UVw);
      }
      if(ray.next < 0) {
        if(prlev) printf("### escape\n");
        break;
      }
      /* if(old >= 0) assert(ray.here == old);  */
      old = ray.next;
      ncross++;
      neps--;
      if(dist >  eps) neps = 10;
      if((prlev > 1) && (neps < 1)) {
        printf("*** %s:%d dist=%.2le iso=%d\n",__FILE__,__LINE__,dist,iso);
        printf(" vac=%d inUV=%d here=%d next=%d", vacuum, geom->InUV(), ray.here,ray.next);
        printf(" ncross=%d",ncross);
        printf(" iter=%d \n",iter);
        printf("  %.3lf %.3lf %.3lf\n",UVp[0],UVp[1],UVp[2]);
      }
      assert(neps >  -10);
    }	// while ncross
    CGp[0] += 0.05;
    CGp[1] += 0.05;

  }	// for iter
  printf("=== done ===\n");
  tally->print();
  delete geom;
  for(int n=0; n < ncgs; n++) delete cg[n];
  delete [] cg;
};
