/*
	cgset.cpp
*/
#include <stdio.h>
#include <math.h>
#include "CGcyl.h"
#include "CGrpp.h"
#include "CGtrc.h"
#include "CGset.h"
#include "Xor128.h"


int main()
{
#ifdef Pro1
  int ncgs = 3;
  CG **cg = new CG*[ncgs];
  cg[0] = new CGcyl();
  cg[0]->CYL(1.0,0.0,10.0);
  cg[1] = new CGrpp();
  cg[1]->RPP(-6.0,6.0, -6.0,6.0, 0.0,10.0);
  cg[2] = new CGrpp();
  cg[2]->RPP(-15.0,15.0, -15.0,15.0, 0.0,10.55);
#else
  int ncgs = 3;
  CG **cg = new CG*[ncgs];
  cg[0] = new CGrpp();
  // cg[0]->RPP(-17.55,17.55, -17.55,17.55, -0.1, 10.56);
  double xmin[] = {-17.55,-17.55,-0.1};
  double xmax[] = { 17.55, 17.55,10.56};
  cg[0]->RPP(xmin,xmax);
  cg[1] = new CGtrc();
  cg[1]->TRC(8.0,17.55,1.0,9.55);
  cg[2] = new CGcyl();
  cg[2]->CYL(8.0, 0.0,1.0);
  // double p[] = {0.0,0.5,0.1};
  // double p[] = {0.0,2.0,5.0};
  // double p[] = {10.0,-10.0, 10.0};
  // double p[] = {0.0,0.0,0.5};
  double p[] = {0.0,0.0, 5.0};
  
#endif
  CGset *cgset = new CGset();
  cgset->set(ncgs, cg);
  double w[] = {0.0,0.0,0.0};
  w[2] = -sqrt(1.0-w[0]*w[0]-w[1]*w[1]);

  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();
  int nran = 100;
  printf("#  CGset\n");
  printf("# dist  px py pz    dx dy dz\n");
  for(int ntry=0; ntry < nran; ntry++) {
    ran->direcS2(w);
    double dist = cgset->track(p,w);
    
    printf("%.4le  %.3lf %.3lf %.3lf  ",dist, p[0],p[1],p[2]);
    printf("%.3lf %.3lf %.3lf\n", w[0]*dist,w[1]*dist,w[2]*dist);

  }	// for nran
  for(int n=0; n < ncgs; n++) {
    delete cg[n];
  }
  delete [] cg;
  delete cgset;
  delete ran;
};