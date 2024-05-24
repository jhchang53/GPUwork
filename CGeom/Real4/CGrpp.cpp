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

void CGrpp::RPP(const float a[], const float b[])
{
  for(int i=0; i < 3; i++) {
    box_C[i] = (b[i]+a[i])/2;
    box_e[i] = (b[i]-a[i])/2;
  }
  if(prlev) {
    printf("C=(%.3lf,%.3lf,%.3lf)",box_C[0],box_C[1],box_C[2]);
    printf(" e=(%.3lf,%.3lf,%.3lf)",box_e[0],box_e[1],box_e[2]);
    printf("\n");
  }
};	// CGrpp::Box

int CGrpp::Clip(float denom, float numer, float &t0, float &t1)
{
  // Liang{Barsky clipping for a linear component against a face of a box.
  if(denom > 0.0) {
    if(numer > denom*t1) return 0;
    if(numer > denom*t0) t0 = numer/denom;
    return 1;
  }
  else if(denom < 0.0) {
    if(numer > denom*t0) return 0;
    if(numer > denom*t1) t1 = numer/denom;
    return 1;
  }
  else {
    if(numer <= 0.0) return 1;
    else return 0;
  }
  return 0;
};      // CGrpp::Clip


int CGrpp::track(const float p[], const float w[], float trk[])
{
  /*  same as Box::DoLineQuery  */
  float P[3];
  for(int i=0; i < 3; i++) P[i] = p[i]-box_C[i];
  float t0 = -INF;
  float t1 = +INF;
  int numPoints;
  int notCulled =
    Clip(+w[0], -P[0]-box_e[0], t0,t1) &&
    Clip(-w[0], +P[0]-box_e[0], t0,t1) &&
    Clip(+w[1], -P[1]-box_e[1], t0,t1) &&
    Clip(-w[1], +P[1]-box_e[1], t0,t1) &&
    Clip(+w[2], -P[2]-box_e[2], t0,t1) &&
    Clip(-w[2], +P[2]-box_e[2], t0,t1);
  if(notCulled) {
    if(t1 > t0) {
      // The intersection is a segment P + t * w with t in [ t0 , t1 ] .
      numPoints = 2;
      trk[0] = t0;
      trk[1] = t1;
    }
    else {
      // The intersection is a segment P + t * w with t = t0 .
      numPoints = 1;
      trk[0] = t0;
      trk[1] = t0;
    }
  }
  else {
    // The line does not intersect the box. Return invalid parameters .
    numPoints = 0;
    trk[0] = +INF;
    trk[1] = -INF;
  }
  if(numPoints < 2) numPoints = 0;	// fix 
  return numPoints;
};	// CGrpp::track
