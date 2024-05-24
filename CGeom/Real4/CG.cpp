/*
	CG.cpp - generic class
*/
#include <stdio.h>
#include <stdlib.h>
#include "CG.h"

CG::CG()
{
  prlev = 0;
  INF = 1.0/0.0;
  eps = 1.0e-20;
};

CG::~CG()
{
};

void CG::CYL(const float R, const float H0, const float H1)
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};

void CG::TRC(const float R0, const float R1, const float Z0, const float Z1)
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};


void CG::RPP(const float bmin[], const float bmax[])
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};

int CG::track(const float p[], const float w[], float trk[])
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};
