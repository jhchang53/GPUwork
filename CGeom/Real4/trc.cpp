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
  float R0 = 1.0;
  float R1 = 2.0;
  float Z0 = 0.5;
  float Z1 = 3.0;
  cgeom->TRC(R0,R1,Z0,Z1);
  float p[] = {0.0,0.0,-1.0};
  printf("# Z0=%.2lf Z1=%.2lf R0=%.2lf R1=%.2lf\n",Z0,Z1,R0,R1);
  printf("# p=(%.3lf,%.3lf,%.3lf)\n",p[0],p[1],p[2]);
  float w[] = {0.316,0.0,1.0};
  w[0] = sqrt(0.1);
  w[2] = sqrt(1.0-w[0]*w[0]-w[1]*w[1]);
  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();
  int nran = 100;
  printf("#  TRC trk  (px,py,pz)  (dx,dy,dz)\n");
  float trk[2];
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
