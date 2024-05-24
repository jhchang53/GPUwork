/*
	box.cpp - find intersection of Aligned box with ray
*/
#include <stdio.h>
#include <stdlib.h>
#include "Box.h"

int main()
{
  Box *box = new Box();
  
  double bmin[] = {0.0,0.0,0.0};
  double bmax[] = {1.0,1.0,2.0};
  box->set(bmin,bmax);
  double Pin[] = {0.0,1.0,-0.1};
  double D[] = {0.0,0.5,0.5};
  printf("Pin: ");
  for(int i=0; i < 3; i++) printf(" %.2lf",Pin[i]);
  printf("   D:");
  for(int i=0; i < 3; i++) printf(" %.2lf",D[i]);
  printf("\n");

  double t[2];
  int np = box->DoRayQuery(Pin,D,t);
  printf("np=%d t=%.2le,%.2le\n",np,t[0],t[1]);
  for(int n=0; n < np; n++) {
    printf("t=%.2lf ",t[n]);
    for(int i=0; i < 3; i++) {
      printf(" %.2lf",Pin[i] + t[n]*D[i]);
    }
    printf("\n");
  }
};
