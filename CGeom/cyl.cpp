/*
	cyl.cpp
*/
#include <stdio.h>
#include <math.h>
#include "CGcyl.h"
#include "Xor128.h"

int main()
{
  CGcyl *cgeom = new CGcyl();
  double R = 1.0;
  double H0 = 0.0;
  double H1 = 2.0;
  cgeom->CYL(R,H0,H1);
  double p[] = {0.0,0.0,0.5};
  double w[] = {0.0,0.0,0.0};

  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();
  double trk[2];
  int nran = 500;
  printf("# CYL tmin,tmax (px,py,pz)  (dx,dy,dz)\n");
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
