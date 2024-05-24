/*
	CGset.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include "CGset.h"

CGset::CGset()
{
  prlev = 0;
  eps = 1.0e-20;
};

CGset::~CGset()
{
};

void CGset::setPrlev(int prl)
{
  prlev = prl;
};	// CGset::setPrlev

void CGset::set(int ncgs_i, CG **cg_i)
{
  ncgs = ncgs_i;
  cg = cg_i;
  isok = new int[ncgs];
  t_min = new float[ncgs];
  t_max = new float[ncgs];
  isin  = new int[ncgs];
  isnext = new int[ncgs];
};

float CGset::track(const float p[], const float w[])
{
  /*  find positive minimum */
  if(prlev) {
    printf("== track p=");
    printf(" %.3lf,%.3lf,%.3lf ",p[0],p[1],p[2]);
    printf(" %.3lf,%.3lf,%.3lf ",w[0],w[1],w[2]);
    printf("\n");
  }
  float tfarmin = 1.0/0.0;
  float trk[2];
  int thiszone = -1;
  for(int n=0; n < ncgs; n++) {
    if(prlev) printf("CG%d:",n);
    CG *cgeom = cg[n];
    int ok = cgeom->track(p,w, trk);
    isok[n] = ok;
    if(!ok) continue;
    float tmin = trk[0];
    float tmax = trk[1];
    t_min[n] = tmin;
    t_max[n] = tmax;
    if(tmax > 0.0) {
      if(tfarmin > tmax) {
        tfarmin = tmax;
        thiszone = n;
      }
    }
    if(prlev) {
      printf(" %.5lf %.5lf ",tmin,tmax);
      float wx = w[0]; float wy = w[1];  float wz = w[2];
      float pxmin = p[0]+tmin*wx;
      float pymin = p[1]+tmin*wy;
      float pzmin = p[2]+tmin*wz;
      printf("  %.3lf %.3lf %.3lf ",pxmin,pymin,pzmin);
      float pxmax = p[0]+tmax*wx;
      float pymax = p[1]+tmax*wy;
      float pzmax = p[2]+tmax*wz;
      printf("  %.3lf %.3lf %.3lf\n",pxmax,pymax,pzmax);
    }
  }   // end of p is  inside * */
  if(thiszone < 0) return 0.0;
  for(int n=0; n < ncgs; n++) {
    if(isok[n]) {
      isin[n] = 0;
      if((t_min[n] <= 0.0) && (t_max[n] >= 0.0)) isin[n] = 1;
      isnext[n] = 0;
      if((t_min[n] <= tfarmin) && (t_max[n] >= tfarmin)) isnext[n] = 1;
    }
    else {
      isin[n] = 0;
      isnext[n] = 0;
    }
  }
  // isin[thiszone] = 1;	// force thiszone as 1
  if(prlev) {
    printf(" thiszone=%d tfarmin=%.5le isin=[",thiszone,tfarmin);
    for(int n=0; n < ncgs; n++) printf(" %d",isin[n]);
    printf("]");
    printf("  next=[");
    for(int n=0; n < ncgs; n++) printf(" %d",isnext[n]);
    printf("]\n");

    for(int n=0; n < ncgs; n++) {
      printf("%d(%d): %.5lf %.5lf\n",n,isok[n],t_min[n],t_max[n]);
    }
  }
  return tfarmin;
};	// CGset::track

int CGset::cgthis()
{
  int zone = 0;
  int p = 1;
  for(int n=0; n < ncgs; n++) {
    if(isin[n]) {
      zone += p;
    }
    p = 2*p;
  }
  return zone;
};	// CGset::cgthis

int CGset::cgnext()
{
  int zone = 0;
  int p = 1;
  for(int n=0; n < ncgs; n++) {
    if(isnext[n]) {
      zone += p;
    }
    p = 2*p;
  }
  return zone;
};      // CGset::cgnext

