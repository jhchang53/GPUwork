/*
   tranrot.cpp - float version
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TRcoord.h"

int main()
{
  TRcoord *tr = new TRcoord();
  float Zb = 20.0;
  float theta = 15.0;
  float phi = 0.0;
  float psi = 30.0;
  float target[] = {5.0,5.0,15.0};
  tr->set(target,Zb,theta,phi,psi);
  float PI = 2*acos(0.0);
  float qvec[3];
  qvec[0] = cos(45.0*PI/180.0);
  qvec[1] = sin(30.0*PI/180.0);
  qvec[2] = sqrt(1.0-qvec[0]*qvec[0]-qvec[1]*qvec[1]);
  float uvec[3];
  tr->S2Uvec(qvec,uvec);
};

