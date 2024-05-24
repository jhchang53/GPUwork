/*
	cgsera.cpp - complicated test
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
/*
  RPP    1      -8.5       8.5    -6.25      6.25      74.0      75.72
  RPP    2     -25.0      25.0    -25.0      25.0      61.5      73.0
  RPP    3     -24.0      24.0    -24.0      24.0      50.0      55.0
  RPP    4     -34.0      34.0    -34.0      34.0      23.5      41.5
  RPP    5     -25.0      25.0    -25.0      25.0      18.5      22.7
  RPP    6     -35.0      35.0    -35.0      35.0      18.5      76.0
  RPP    7     -35.1      35.1    -35.1      35.1       0.0      18.5
  END
    1    1     +1
    2    2     +2
    3    3     +3
    4    4     +4
    5    5     +5
    6    6     +6     -1     -2     -3     -4     -5
    7    7     +7     -6
  END
*/
  int ncgs = 7;
  CG **cg = new CG*[ncgs];
  cg[0] = new CGrpp();  	// we set RPP6 as base
  // cg[0]->RPP( -35.0,     35.0,   -35.0,     35.0,     18.5,     76.0);
  double xmin0[] = {-35.0, -35.0, 18.5};
  double xmax0[] = { 35.0,  35.0, 76.0};
  cg[0]->RPP(xmin0,xmax0);
  cg[1] = new CGrpp();
  // cg[1]->RPP( -8.5 ,     8.5 ,   -6.25,     6.25,     74.0,     75.72);
  double xmin1[] = { -8.5 ,  -6.25, 74.0};
  double xmax1[] = {  8.5 ,   6.25, 75.72};
  cg[1]->RPP(xmin1,xmax1);
  cg[2] = new CGrpp();
  // cg[2]->RPP( -25.0,     25.0,   -25.0,     25.0,     61.5,     73.0);
  double xmin2[] = {-25.0, -25.0, 61.5};
  double xmax2[] = { 25.0,  25.0, 73.0};
  cg[2]->RPP(xmin2,xmax2);
  cg[3] = new CGrpp();
  // cg[3]->RPP( -24.0,     24.0,   -24.0,     24.0,     50.0,     55.0);
  double xmin3[] = {-24.0, -24.0, 50.0};
  double xmax3[] = { 24.0,  24.0, 55.0};
  cg[3]->RPP(xmin3,xmax3);
  cg[4] = new CGrpp();
  // cg[4]->RPP( -34.0,     34.0,   -34.0,     34.0,     23.5,     41.5);
  double xmin4[] = {-34.0, -34.0, 23.5};
  double xmax4[] = { 34.0,  34.0, 41.5};
  cg[4]->RPP(xmin4,xmax4);
  cg[5] = new CGrpp();
  // cg[5]->RPP( -25.0,     25.0,   -25.0,     25.0,     18.5,     22.7);
  double xmin5[] = {-25.0, -25.0, 18.5};
  double xmax5[] = { 25.0,  25.0, 22.7};
  cg[5]->RPP(xmin5,xmax5);
  cg[6] = new CGrpp();
  // cg[6]->RPP( -35.1,     35.1,   -35.1,     35.1,      0.0,     18.5);
  double xmin6[] = {-35.1, -35.1, 0.0};
  double xmax6[] = { 35.1,  35.1, 18.5};
  cg[6]->RPP(xmin6,xmax6);
  int nbound = 2;
  int bounding[] = {1,64};
  int nzone = 1;
  for(int j=0; j < ncgs; j++) nzone = 2*nzone;
  /*  assign material  */
  int *reg = new int[nzone];
  for(int i=0; i < nzone; i++) reg[i] = -1;
  /*        m      0123456	*/
  reg[1]  = 0;	// 1000000
  reg[3]  = 1;	// 1100000
  reg[5]  = 2;	// 1010000
  reg[9]  = 3;	// 1001000
  reg[17] = 4;	// 1000100
  reg[33] = 5;	// 1000010
  reg[64] = 6;	// 0000001
  reg[65] = 7;	// 1000001
  reg[97] = 8;	// 1000011

  double p0[] = {0.0,0.0,75.9};

  CGset *cgset = new CGset();
  cgset->set(ncgs, cg);
  // cgset->setPrlev(1);
  double w[] = {0.0,0.0,1.0};

  /*  random number */
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 20; nj++) ran->jump();
  int nran = 100;
  printf("#  CGsera\n");
  printf("# i m j dist  px py pz    dx dy dz\n");
  for(int ntry=0; ntry < nran; ntry++) {
    printf(" ntry=%d\n",ntry);

    // ran->direcS2(w);
    double p[3];
    for(int j=0; j < 3; j++) p[j] = p0[j];
    for(int ngo=0; ngo < 1000; ngo++) {
      ran->direcS2(w);
      double dist = cgset->track(p,w);
      int cgthis = cgset->cgthis();
      int cgnext = cgset->cgnext();
      int mat = reg[cgthis];
      if(mat < 0) {
        printf("*** cgthis=%d\n",cgthis);
        exit(0);
      }
      int escape = 0;
      if(dist < 1.0e-20) {
        for(int nb=0;!escape && (nb < nbound); nb++) {
          if(cgnext == bounding[nb]) escape = 1;
        }
      }
      if(escape) break;
      // if(cgnext >= 0) {
        printf("%d %d %d ",cgthis,mat,cgnext);
        printf("%.4le  %.3lf %.3lf %.3lf  ",dist, p[0],p[1],p[2]);
        printf("%.3lf %.3lf %.3lf", w[0]*dist,w[1]*dist,w[2]*dist);
        p[0] += w[0]*dist;   p[1] += w[1]*dist;  p[2] += w[2]*dist;
        // printf(" => (%.3lf) %.3lf %.3lf %.3lf",dist,p[0],p[1],p[2]);
        printf("\n");
      // }
    }
  }	// for nran
  delete [] reg;
  delete cgset;
  for(int n=0; n < ncgs; n++) {
    delete cg[n];
  }
  delete [] cg;
  delete ran;
};
