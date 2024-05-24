/*
	vox.cpp	- float version
*/
#include <stdio.h>
#include <math.h>
#include "UVfile.h"
#include "Box.h"
#include "Vox.h"

int main()
{
  int nxyz[] = {128,128,134};
  float dxyz[] = {0.2,0.2,0.2};
  UVfile *uvf = new UVfile();
  uvf->open("/home/jhchang/Dicom3/test3/0831150844478.uv",nxyz,dxyz);
  float bmin[] = {0.0,0.0,0.0};
  float bmax[3];
  for(int i=0; i < 3; i++) {
    bmax[i] = nxyz[i]*dxyz[i];
  }
  Box *box = new Box();
  box->set(bmin,bmax);
  Vox *vox = new Vox();
  vox->setLattice(nxyz,dxyz,uvf->UV());
  /*  start  */
  // float P[] = {11.95,12.12,12.14};
  float P[] = {10.0,13.5,-1.0};
  float wxyz[] = {-0.1,0.1,0.0};
  wxyz[2] = sqrt(1.0-wxyz[0]*wxyz[0]-wxyz[1]*wxyz[1]);
  float t[2];
  int np = box->DoRayQuery(P,wxyz,t);
  printf("np=%d t=%.2le,%.2le\n",np,t[0],t[1]);
  for(int n=0; n < np; n++) {
    printf("t=%.2lf ",t[n]);
    for(int i=0; i < 3; i++) {
      printf(" %.2lf",P[i] + t[n]*wxyz[i]);
    }
    printf("\n");
  }
  /*  entrance point */
  printf("#  P=(%.3lf,%.3lf,%.3lf\n",P[0],P[1],P[2]);

  float P0[3];
  for(int i=0; i < 3; i++) P0[i] = P[i] + t[0]*wxyz[i];
  char chx = ' ';
  if(np == 2) {
    for(int n=0; n < 9; n++) {
      Ray ray = vox->track(P0,wxyz, chx);
      // printf(" chx=%c chn=%c track=%.3lf\n",ray.chx,ray.chn,ray.dist);
      float dist = ray.dist;
      if((ray.next < 1) || (dist < 1.0e-20)) break;
      for(int i=0; i < 3; i++) P0[i] += dist*wxyz[i];
      chx = ray.next;
    }
  };
};
