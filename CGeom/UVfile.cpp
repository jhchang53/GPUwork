/*
	UVfile.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include "UVfile.h"

UVfile::UVfile()
{
  uv = nullptr;
};

UVfile::~UVfile()
{
  delete [] uv;
};

void UVfile::open(const char *uvpath, const int nxyz[], const double dxyz[])
{
  FILE *F = fopen(uvpath,"r");
  if(F == NULL) {
    printf("*** %s not exist.\n",uvpath);
    exit(0);
  }
  nx = nxyz[0];	ny = nxyz[1];	nz = nxyz[2];
  long nuvs = nx*ny*nz;
  uv = new char[nuvs];
  fgets(uv,nuvs,F);
  fclose(F);
  int pbeg = 0;
  while(uv[pbeg] == ' ') pbeg++;
  pbeg += 10*nx*ny + 3*nx;
};	// UVfile::open
