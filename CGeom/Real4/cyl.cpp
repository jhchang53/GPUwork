/*
	cyl.cpp	- flaot version
*/
#include <stdio.h>
#include <math.h>
#include "CGcyl.h"
#include "Xor128.h"

int main()
{
  CGcyl *cgeom = new CGcyl();
  float R = 1.0;
  float H0 = 0.0;
  float H1 = 2.0;
  cgeom->CYL(R,H0,H1);
  float p[] = {0.0,0.0,0.5};
  float w[] = {0.0,0.0,0.0};

  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();
  float trk[2];
  int nran = 500;
  printf("# CYL tmin,tmax (px,py,pz)  (dx,dy,dz)\n");
  for(int n=0; n < nran; n++) {
    ran->direcS2(w);
    int ok = cgeom->track(p,w, trk);
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
      printf("  %.3lf %.3lf %.3lf\n",pxmax-pxmin,pymax-pymin,pzmax-pzmin);
    }
  }
};
