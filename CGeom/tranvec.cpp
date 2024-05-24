/*
   tranrot.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TRcoord.h"

int main()
{
  TRcoord *tr = new TRcoord();
  double Zb = 20.0;
  double theta = 15.0;
  double phi = 0.0;
  double psi = 30.0;
  double target[] = {5.0,5.0,15.0};
  tr->set(target,Zb,theta,phi,psi);
  double PI = 2*acos(0.0);
  double qvec[3];
  qvec[0] = cos(45.0*PI/180.0);
  qvec[1] = sin(30.0*PI/180.0);
  qvec[2] = sqrt(1.0-qvec[0]*qvec[0]-qvec[1]*qvec[1]);
  double uvec[3];
  tr->S2Uvec(qvec,uvec);
  delete tr;
};

