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

void CG::CYL(const double R, const double H0, const double H1)
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};

void CG::TRC(const double R0, const double R1, const double Z0, const double Z1)
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};


void CG::RPP(const double bmin[], const double bmax[])
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};

int CG::track(const double p[], const double w[], double trk[])
{
  printf("*** %s:%d \n",__FILE__,__LINE__);
  exit(0);
};
