/*
	elas.cpp - isotropic elastic scattering
*/
#include <stdio.h>
#include <math.h>
#include "Xor128.h"
#include "Scatter.h"

int main()
{
  init_xor128();

  Scatter *scat = new Scatter();
  // double w[] = {0.5,0.5,0.0};
  // double w[] = {0.6,0.1,0.0};
  // scat->setPrlev(1);
  double w[] = {0.6,0.5,0.0};
  w[2] = sqrt(1.0-w[0]*w[0]-w[1]*w[1]);
  double A = 12.0;
  double E = 1.0;
  for(int n=0; n < 5; n++) {
    double En = scat->isotropic(A,E, w);
    printf(" E=%.4le En=%.4le\n",E,En);
    E = En;
  }
  delete scat;
};
