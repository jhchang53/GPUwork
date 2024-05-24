/*
	trc.cpp
*/
#include <stdio.h>
#include <math.h>
#include "CGtrc.h"
#include "Xor128.h"

int main()
{
  CGtrc *cgeom = new CGtrc();
  double R0 = 1.0;
  double R1 = 2.0;
  double Z0 = 0.5;
  double Z1 = 3.0;
  cgeom->TRC(R0,R1,Z0,Z1);
  double p[] = {0.0,0.0,-1.0};
  printf("# Z0=%.2lf Z1=%.2lf R0=%.2lf R1=%.2lf\n",Z0,Z1,R0,R1);
  printf("# p=(%.3lf,%.3lf,%.3lf)\n",p[0],p[1],p[2]);
  double w[] = {0.316,0.0,1.0};
  w[0] = sqrt(0.1);
  w[2] = sqrt(1.0-w[0]*w[0]-w[1]*w[1]);
  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();
  int nran = 100;
  printf("#  TRC trk  (px,py,pz)  (dx,dy,dz)\n");
  double trk[2];
  for(int n=0; n < nran; n++) {
    ran->direcS2(w);
    int ok = cgeom->track(p,w, trk);
    if(ok) {
      double tmin = trk[0];
      double tmax = trk[1];
      printf(" %.4lf %.4lf ",tmin,tmax);
      double wx = w[0]; double wy = w[1];  double wz = w[2];
      double pxmin = p[0]+tmin*wx;
      double pymin = p[1]+tmin*wy;
      double pzmin = p[2]+tmin*wz;
      printf("  %.3lf %.3lf %.3lf ",pxmin,pymin,pzmin);
      double pxmax = p[0]+tmax*wx;
      double pymax = p[1]+tmax*wy;
      double pzmax = p[2]+tmax*wz;
      printf("  %.3lf %.3lf %.3lf\n",pxmax-pxmin,pymax-pymin,pzmax-pzmin);
    }
  }
  delete ran;
  delete cgeom;
};
