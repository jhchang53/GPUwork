/*
 	BOX.cu
	ref) D.Eberly, Intersection of a Line and a Box, Geometric Tools, 2018
*/
#include <stdio.h>

__device__ double box_C[3],box_e[3],box_min[3],box_max[3];

__device__ void BOX_set(double a[], double b[])
{
  for(int i=0; i < 3; i++) {
    box_C[i] = (b[i]+a[i])/2;
    box_e[i] = (b[i]-a[i])/2;
    box_min[i] = a[i];
    box_max[i] = b[i];
  }
};      // BOX_set

__device__ int BOX_isInside(double P[])
{
  /*  check whether a point is inside of the box  */
  if(P[0] < box_min[0]) return 0;
  if(P[0] > box_max[0]) return 0;
  if(P[1] < box_min[1]) return 0;
  if(P[1] > box_max[1]) return 0;
  if(P[2] < box_min[2]) return 0;
  if(P[2] > box_max[2]) return 0;
  return 1;
};      // BOX_isInside

__device__ int BOX_Clip(double denom, double numer, double &t0, double &t1)
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
};      // BOX_Clip

__device__ int BOX_DoLineQuery(const double P[], const double D[], double t[2])
{
  /*  DoLineQuery is called with size 3 arrays of Position and Directtion */
  double INF = 1.0/0.0;
  double t0 = -INF;
  double t1 = +INF;
  int numPoints;
  int notCulled =
    BOX_Clip(+D[0], -P[0]-box_e[0], t0,t1) &&
    BOX_Clip(-D[0], +P[0]-box_e[0], t0,t1) &&
    BOX_Clip(+D[1], -P[1]-box_e[1], t0,t1) &&
    BOX_Clip(-D[1], +P[1]-box_e[1], t0,t1) &&
    BOX_Clip(+D[2], -P[2]-box_e[2], t0,t1) &&
    BOX_Clip(-D[2], +P[2]-box_e[2], t0,t1);
  if(notCulled) {
    if(t1 > t0) {
      // The intersection is a segment P + t * D with t in [ t0 , t1 ] .
      numPoints = 2;
      t[0] = t0;
      t[1] = t1;
    }
    else {
      // The intersection is a segment P + t * D with t = t0 .
      numPoints = 1;
      t[0] = t0;
      t[1] = t0;
    }
  }
  else {
    // The line does not intersect the box. Return invalid parameters .
    numPoints = 0;
    t[0] = +INF;
    t[1] = -INF;
  }
  return numPoints;
};      // Box::DoLineQuery


__device__ int BOX_DoRayQuery(const double Pin[], double t[2])
{
  /*  DoRayQuery is called with size 8 array  */
  double p[3],w[3];
  for(int i=0; i < 3; i++) {
   p[i] = Pin[2+i]-box_C[i];
   w[i] = Pin[5+i];
  }
  int numPoints = BOX_DoLineQuery(p,w,t);
  if(numPoints > 0) {
    // The line containing the ray intersects the box in the interval
    // [t0,t1] . Compute the intersection of [t0,t1] and [0,+infinity).
    if(t[1] >= 0) {
      if(t[0] < 0) t[0] = 0.0;
    }
    else numPoints = 0;
  }
  return numPoints;
};      // Box::DoRayQuery

