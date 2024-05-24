/*
	rpp.cpp
*/
#include <stdio.h>
#include <math.h>
#include "CGrpp.h"
#include "Xor128.h"

int main()
{
  CGrpp *cgeom = new CGrpp();
  double xmin[] = {1.0,1.5,2.0};
  double xmax[] = {2.0,2.5,3.0};
  // cgeom->RPP(xmin[0],xmax[0], xmin[1],xmax[1], xmin[2],xmax[2] );
  cgeom->RPP(xmin,xmax);
  double p[] = {1.5,2.0,1.8};
  double w[] = {0.5,0.5,0.0};
  w[2] = sqrt(1.0-w[0]*w[0]-w[1]*w[1]);
  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();
  double trk[2];
  printf("#  tmin tmax  (px,py,pz)  (dx,dy,dz)\n");
  int nran = 1000;
  for(int n=0; n < nran; n++) {
    ran->direcS2(w);
    // printf("  w=%.6lf\n",w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    int ok = cgeom->track(p, w, trk);
    // printf(" ok==%d\n",ok);
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
      printf("  %.3lf %.3lf %.3lf \n",pxmax-pxmin,pymax-pymin,pzmax-pzmin);
    }
  }
  delete ran;
  delete cgeom;
};
