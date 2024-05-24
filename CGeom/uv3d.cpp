/*
	uv3d.cpp
*/
#include <stdio.h>
#include <math.h>
#include "UVfile.h"
#include "Vox.h"

int main()
{
  int nxyz[] = {128,128,134};
  double dxyz[] = {0.2,0.2,0.2};
  UVfile *uvf = new UVfile();
  uvf->open("/home/jhchang/Dicom3/test3/0831150844478.uv",nxyz,dxyz);
  Vox *vox = new Vox();
  double xyz0[] = {0.0,0.0,0.0};
  vox->setLattice(nxyz,dxyz,uvf->UV());
  double pxyz[] = {11.95,12.12,12.14};
  double wxyz[] = {-0.5,-0.5,0.0};
  wxyz[2] = -sqrt(1.0-wxyz[0]*wxyz[0]-wxyz[1]*wxyz[1]);
  char chx = '9';
  Ray ray = vox->track(pxyz,wxyz,chx);
  printf(" ray here=%c next=%c dist=%.3lf\n",ray.here,ray.next,ray.dist);
  delete vox;
  delete uvf;
};
