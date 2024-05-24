/*
	inel.cpp - inelastic scattering
*/
#include <stdio.h>
#include <math.h>
#include "Xor128.h"
#include "Scatter.h"

int main()
{
  init_xor128();

  Scatter *scat = new Scatter();
  scat->setPrlev(0);
  double w[] = {0.6,0.5,0.0};
  w[2] = sqrt(1.0-w[0]*w[0]-w[1]*w[1]);
  double A = 12.0;
  double p0 = 1.0;  double p1 = 0.0;
  double wvec[3];
  int bin[100];
  for(int j=0; j < 100; j++) bin[j] = 0.0;
  for(int n=0; n < 100000; n++) {
    double xmu = scat->P1(p0,p1, w, wvec);
    int b = 100*(xmu+1)/2;
    bin[b%100]++;
  }
  delete scat;
  for(int i=0; i < 100; i++) printf("%d\n",bin[i]);
};
