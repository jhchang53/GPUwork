/*
	rpp.cpp - float version
*/
#include <stdio.h>
#include <math.h>
#include "CGrpp.h"
#include "Xor128.h"

int main()
{
  CGrpp *cgeom = new CGrpp();
  float xmin[] = {1.0,1.5,2.0};
  float xmax[] = {2.0,2.5,3.0};
  cgeom->RPP(xmin,xmax);
  float p[] = {1.5,2.0,1.8};
  float w[] = {0.5,0.5,0.0};
  w[2] = sqrt(1.0-w[0]*w[0]-w[1]*w[1]);
  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();
  float trk[2];
  printf("#  tmin tmax  (px,py,pz)  (dx,dy,dz)\n");
  int nran = 1000;
  for(int n=0; n < nran; n++) {
    ran->direcS2(w);
    // printf("  w=%.6lf\n",w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    int ok = cgeom->track(p, w, trk);
    // printf(" ok==%d\n",ok);
    if(ok) {
      float tmin = trk[0];
      float tmax = trk[1];
      printf(" %.4lf %.4lf ",tmin,tmax);
      float wx = w[0]; float wy = w[1];  float wz = w[2];
      float pxmin = p[0]+tmin*wx;
      float pymin = p[1]+tmin*wy;
      float pzmin = p[2]+tmin*wz;
      printf("  %.3lf %.3lf %.3lf ",pxmin,pymin,pzmin);
      float pxmax = p[0]+tmax*wx;
      float pymax = p[1]+tmax*wy;
      float pzmax = p[2]+tmax*wz;
      printf("  %.3lf %.3lf %.3lf \n",pxmax-pxmin,pymax-pymin,pzmax-pzmin);
    }
  }
};
