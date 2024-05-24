/*
	compton.cpp - Compton scattering
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
  double En = 5.0e+6;	// eV
  for(int n=0; n < 1; n++) {
    double xmu = scat->Compton(En);
  }
  delete scat;
};
