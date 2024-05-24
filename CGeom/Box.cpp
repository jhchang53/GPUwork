/*
	Box.cpp
	ref) D.Eberly, Intersection of a Line and a Box, Geometric Tools, 2018
*/
#include <stdio.h>
#include "Box.h"

Box::Box()
{
  prlev = 0;
  INF = 1.0/0.0;
};

Box::~Box()
{
};

void Box::set(double a[], double b[])
{
  for(int i=0; i < 3; i++) {
    box_C[i] = (b[i]+a[i])/2;
    box_e[i] = (b[i]-a[i])/2;
    xmin[i] = a[i];
    xmax[i] = b[i];
  }
};	// Box::set

int Box::isInside(double P[])
{
  /*  check whether a point is inside of the box  */
  if(P[0] < xmin[0]) return 0;
  if(P[0] > xmax[0]) return 0;
  if(P[1] < xmin[1]) return 0;
  if(P[1] > xmax[1]) return 0;
  if(P[2] < xmin[2]) return 0;
  if(P[2] > xmax[2]) return 0;
  return 1;
};	// Box::isInside

int Box::Clip(double denom, double numer, double &t0, double &t1)
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
};	// Box::Clip

int Box::DoLineQuery(const double P[], const double D[], double t[2])
{
  double t0 = -INF;
  double t1 = +INF;
  int numPoints;
  int notCulled = 
    Clip(+D[0], -P[0]-box_e[0], t0,t1) &&
    Clip(-D[0], +P[0]-box_e[0], t0,t1) &&
    Clip(+D[1], -P[1]-box_e[1], t0,t1) &&
    Clip(-D[1], +P[1]-box_e[1], t0,t1) &&
    Clip(+D[2], -P[2]-box_e[2], t0,t1) &&
    Clip(-D[2], +P[2]-box_e[2], t0,t1);
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
};	// Box::DoLineQuery


int Box::DoRayQuery(const double Pin[], const double D[], double t[2])
{
  double P[3];
  for(int i=0; i < 3; i++) P[i] = Pin[i]-box_C[i];
  int numPoints = DoLineQuery(P,D,t);
  if(numPoints > 0) {
    // The line containing the ray intersects the box in the interval
    // [t0,t1] . Compute the intersection of [t0,t1] and [0,+infinity).
    if(t[1] >= 0) {
      if(t[0] < 0) t[0] = 0.0;
    }
    else numPoints = 0;
  }
  return numPoints;
};	// Box::DoRayQuery
