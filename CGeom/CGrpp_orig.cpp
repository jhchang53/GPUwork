/*
	CGrpp.cpp
	ref) P.J.Schneider and D.H. Eberly, Geometric Tools for Computer Graphics,
	 Elsevier, 2003. 
*/
#include <stdio.h>
#include <limits>
#include <math.h>
#include <assert.h>
#include "CGrpp.h"

using namespace std;

CGrpp::CGrpp()
{
};

CGrpp::~CGrpp()
{
};

void CGrpp::RPP(double x0, double x1, double y0, double  y1, double z0, double z1)
{
  xmin[0] = x0;   xmin[1] = y0;  xmin[2] = z0;
  xmax[0] = x1;   xmax[1] = y1;  xmax[2] = z1;
  for(int n=0; n < 3; n++) {
    assert(xmax[n] > xmin[n]);
  }
};	// CGrpp::RPP

int CGrpp::track(double p[], double w[], double trk[])
{
  /*  Schneider chap.11.12.2  */
  /*  Linear Component and Axis-Aligned Bounding Box (AABB)	*/
  /*  Schneider, Chap.11.12.2	*/
  int prlev = 1;
  if(prlev) {
    printf("*** RPP  p=(%.3lf %.3lf %.3lf)",p[0],p[1],p[2]);
    printf("  w=(%.3lf %.3lf %.3lf)\n",w[0],w[1],w[2]);
  }
  double tmax = INF;
  double tmin = -INF;
  for(int n=0; n < 3; n++) {	// for each pair of plnaes
    assert(xmax[n] > xmin[n]);
    if(prlev) printf("  x=[%.3lf  %.3lf]  wx=%.4lf\n",xmin[n],xmax[n], w[n]);
    if(fabs(w[n]) < eps) {
      // ray parallel to planes
      // if((p[n] < xmin[n]) || (p[n] < xmax[n])) return 0;
      continue;
    }
    double t0 = (xmin[n]-p[n])/w[n];
    double t1 = (xmax[n]-p[n])/w[n];
    if(prlev) printf("oo t0=%.3lf t1=%.3lf\n",t0,t1);

    // make t0 smaller
    if(t0 > t1) {	// swap
      double tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    if(prlev) {
      printf("t0=%.3lf (%.3lf,%.3lf,%.3lf) ",t0,p[0]+t0*w[0],p[1]+t0*w[1],p[2]+t0*w[2]);
      printf("t1=%.3lf (%.3lf,%.3lf,%.3lf) ",t1,p[0]+t1*w[0],p[1]+t1*w[1],p[2]+t1*w[2]);
      printf("\n");
    }
    // compare
    if(t0 > tmin) tmin = t0;
    if(t1 < tmax) tmax = t1;
    if(prlev) printf(" tmin=%.3lf tmax=%.3lf\n",tmin,tmax);
  }
  trk[0] = tmin;
  trk[1] = tmax;
  if(tmin > tmax) return 0;

  double x0 = p[0] + tmin*w[0];
  double y0 = p[1] + tmin*w[1];
  double z0 = p[2] + tmin*w[2];
  if(prlev) printf(" P0=(%.3lf,%.3lf,%.3lf) ",x0,y0,z0);
  double x1 = p[0] + tmax*w[0];
  double y1 = p[1] + tmax*w[1];
  double z1 = p[2] + tmax*w[2];
  if(prlev) printf(" P1=(%.3lf,%.3lf,%.3lf)\n",x1,y1,z1);
  return 1;
};	// CGrpp::RPP
