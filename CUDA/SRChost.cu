/*
	SRChost.cu
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double *SRCen,*SRCcdf;

int SRCinit_host()
{
  /*  initialize source  */
  double En[] = { 	//  neutron enegrgy (eV)
  1.00000e+07, 7.78800e+06, 6.06530e+06, 4.72370e+06, 3.67880e+06, 2.86500e+06
, 2.23130e+06, 1.73770e+06, 1.35340e+06, 1.05400e+06, 8.20850e+05, 6.39280e+05
, 4.97870e+05, 3.87740e+05, 3.01970e+05, 2.35180e+05, 1.83160e+05, 1.42640e+05
, 1.11090e+05, 8.65200e+04, 6.73800e+04, 5.24700e+04, 4.08700e+04, 3.18300e+04
, 2.47900e+04, 1.93000e+04, 1.50300e+04, 1.17100e+04, 9.12000e+03, 7.10000e+03
, 5.53000e+03, 4.31000e+03, 3.35000e+03, 2.61000e+03, 2.03000e+03, 1.58000e+03
, 1.23000e+03, 9.61120e+02, 7.48520e+02, 6.14420e+01};
  double  spect[] = {	// neutron spectrum
  3.50283e-06, 9.49937e-05, 2.11839e-04, 3.35214e-04, 4.07159e-04, 4.06490e-04
, 3.77376e-04, 3.36183e-04, 3.07504e-04, 2.54209e-04, 1.86235e-04, 1.34501e-04
, 9.97904e-05, 6.42947e-05, 4.11489e-05, 2.93050e-05, 2.20885e-05, 1.64738e-05
, 1.22263e-05, 9.12809e-06, 6.88864e-06, 5.14428e-06, 3.83756e-06, 2.72307e-06
, 2.12038e-06, 1.39512e-06, 1.04543e-06, 7.60837e-07, 4.73659e-07, 3.44188e-07
, 2.34793e-07, 1.68728e-07, 1.19882e-07, 8.29259e-08, 6.09660e-08, 3.09275e-08
, 2.59211e-08, 2.49288e-08, 1.70000e-08};
  int nsrcgrp = sizeof(spect)/sizeof(double);
  assert(sizeof(En)/sizeof(double) == nsrcgrp+1);
  /*  make CDF  */
  double sumE = 0.0;
  for(int n=0; n < nsrcgrp; n++) {
    sumE += spect[n]*(En[n]-En[n+1]);
  };
  SRCen = new double[nsrcgrp+1];
  SRCcdf = new double[nsrcgrp];
  double cdf0 = 0.0;
  for(int n=0; n < nsrcgrp; n++) {
    SRCen[n] = En[n];
    double sE = spect[n]*(En[n]-En[n+1]);
    SRCcdf[n] = cdf0 + sE/sumE;
    cdf0 = SRCcdf[n];
  }
  SRCen[nsrcgrp] = En[nsrcgrp];
  for(int n=0; n < nsrcgrp; n++) {
    printf(" %.5lf",SRCcdf[n]);
  }
  printf("\n");
  printf(" last=%.2le\n",cdf0-1.0);
  return nsrcgrp;
};
