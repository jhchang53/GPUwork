/*
	geom.cpp	- test driver for float version
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "CGcyl.h"
#include "CGrpp.h"
#include "CGtrc.h"
#include "CGset.h"
#include "Ray.h"
#include "Geometry.h"
#include "Xor128.h"

int main()
{
  int ncgs = 3;
  CG **cg = new CG*[ncgs];
  cg[0] = new CGrpp();
  float bmin[] = {-17.55,-17.55,-0.1};
  float bmax[] = { 17.55, 17.55,10.56};
  cg[0]->RPP(bmin,bmax);
  cg[1] = new CGtrc();
  cg[1]->TRC(8.0,17.55,1.0,9.55);
  cg[2] = new CGcyl();
  cg[2]->CYL(8.0, 0.0,1.0);

  Geometry *geom = new Geometry();
  geom->setPrlev(2);
  geom->setCG(ncgs,cg);
  float CGp[] = {7.99,0.0,-1.0};
  float CGw[] = {0.1,0.1,1.0};
  CGw[2] = sqrt(1.0-CGw[0]*CGw[0]-CGw[1]*CGw[1]);
  int old = -1;
  for(int n=0; n < 10; n++) {
    printf(" CGp=(%.3f,%.3f,%.3f)\n",CGp[0],CGp[1],CGp[2]);
    Ray ray = geom->CGtrack(CGp,CGw);
    float dist = ray.dist;
    for(int j=0; j < 3; j++) CGp[j] += dist*CGw[j];
    if(ray.next == 0) break;
    // if(old >= 0) assert(ray.here == old);
      
    old = ray.next;
  }
  printf(" CGp=(%.3f,%.3f,%.3f)\n",CGp[0],CGp[1],CGp[2]);
  delete geom;
  for(int n=0; n < ncgs; n++) delete cg[n];
  delete [] cg;
};
