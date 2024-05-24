/*
	collide.cpp
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
#include "Xor128.h"

int main()
{
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
  double dxyz[] = {0.2,0.2,0.2};
  UVfile *uvf = new UVfile();
  uvf->open("/home/jhchang/Dicom3/test3/0831150844478.uv",nxyz,dxyz);
  geom->setLattice(nxyz,dxyz,uvf->UV());

  TRcoord *trc = new TRcoord();
  double Zb = 20.0; // 20.0;
  double theta = 15.0;
  double phi = 0.0;
  double psi = 0.0;
  double target[] = {12.0,12.0,25.0};
  trc->set(target,Zb,theta,phi,psi);

  /*  set source coordinate  */
  double CGp[] = {0.0,0.0,-1.0};
  double CGw[] = {0.00,0.0,1.0};
  CGw[2] = sqrt(1.0-CGw[0]*CGw[0]-CGw[1]*CGw[1]);

  double UVp[3],UVw[3];
  CGp[0] = -8.0;   CGp[1] = CGp[0];

  init_xor128();
  for(int nj=0; nj < 3; nj++) jump_xor128();

  for(int iter=0; iter < 5; iter++) {
    CGp[0] += 0.5;
    CGp[1] += 0.5;
    printf("### CGp=(%.3lf,%.3lf,%.3lf) ",CGp[0],CGp[1],CGp[2]);
    printf("    CGw=(%.3lf,%.3lf,%.3lf) ",CGw[0],CGw[1],CGw[2]);
    printf("\n");
    trc->S2U(CGp,UVp);
    trc->S2Uvec(CGw,UVw);
    geom->setTR(trc);
    int old = -1;
    double eps = 1.0e-10;
    int ncross = 0;

    while(1) {
      // printf("=== UVp=(%.3lf,%.3lf,%.3lf) ",UVp[0],UVp[1],UVp[2]);
      // printf("   UVw=(%.3lf,%.3lf,%.3lf)\n",UVw[0],UVw[1],UVw[2]);
      Ray ray = geom->track(UVp,UVw);
      if(ray.here < 0) {
        printf("### escape\n");
        break;
      }
      double dist = ray.dist;
      double fflight = 1.0e+20;
      int vacuum = 0;
      if(ray.here == 0) vacuum = 1;
      else if(ray.here == ' ') vacuum = 1;
      if(!vacuum) {
        /*  compute flight distance  -ln(xsi)/Sig_t  */
        double sigtot = 1.0;
        fflight = -log(ranf())/sigtot;
      }
      if(fflight < dist) {
        dist = fflight;
        ray.next = ray.here;
        printf("# collide dist=%.2le\n",dist);
      }
      printf("%d %.3lf %.3lf %.3lf  %.3le",ray.here,UVp[0],UVp[1],UVp[2], dist);
      printf("  %.3lf %.3lf %.3lf",dist*UVw[0],dist*UVw[1],dist*UVw[2]);
      // printf("  here=%d next=%d",ray.here,ray.next);
      printf("  %.3lf %.3lf %.3lf",UVp[0]+dist*UVw[0],UVp[1]+dist*UVw[1],UVp[2]+dist*UVw[2]);
      printf("\n");
      if(dist <= eps) dist = eps;
      for(int j=0; j < 3; j++) UVp[j] += dist*UVw[j];
      // printf("%d ",ray.here);
      // for(int j=0; j < 3; j++) printf(" %.3lf",UVp[j]);
      // printf("\n");
      if(ray.next < 0) {
        printf("### escape\n");
        break;
      }
      // printf("#  old=%d here=%d next=%d dist=%.3lf\n",old,ray.here,ray.next,ray.dist);
      /* if(old >= 0) assert(ray.here == old);  */
      old = ray.next;
      ncross++;
      assert(ncross < 100);
    }	// while ncross
  }	// for iter
  delete geom;
  for(int n=0; n < ncgs; n++) delete cg[n];
  delete [] cg;
  delete trc;
  delete uvf;
};
