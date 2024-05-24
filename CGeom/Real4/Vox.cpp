/*
	Vox.cpp - float version
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "Vox.h"

Vox::Vox()
{
  prlev = 0;
  INF = 1.0e+20;
  eps = 1.0e-20;
  // X = Y = Z = -1;
};

Vox::~Vox()
{
};

void Vox::setLattice(const int nxyz[], const float dxyz[], char* uv_i)
{
  nx = nxyz[0];  dx = dxyz[0];
  ny = nxyz[1];  dy = dxyz[1];
  nz = nxyz[2];  dz = dxyz[2];
  printf(" Box: %.3lf,%.3lf,%.3lf\n",nx*dx,ny*dy,nz*dz);
  
  uv = uv_i;
};	// Vox::setLattice

char Vox::getChar(const float pxyz[])
{
  /*  retrieve char at pxyz, if outside returns 0  */
  int X = pxyz[0]/dx;
  if((X < 0) || (X >= nx)) return 0;
  int Y = pxyz[1]/dy;
  if((Y < 0) || (Y >= ny)) return 0;
  int Z = pxyz[2]/dz;
  if((Z < 0) || (Z >= nz)) return 0;
  return uv[(Z*ny+Y)*nx+X];
};	// Vox::getChar

/*		*/
/*  returns distance of cells with same mat (chx)	*/
/*  until other mat is meet or outside  		*/
/*  when no intersection  returns dist as INF		*/
Ray Vox::track(const float pxyz[], const float wxyz[], char chx)
{
  /* orgin coordinate of lattice is (0,0,0)  */
  float ux = pxyz[0]; int X = ux/dx;
  float uy = pxyz[1]; int Y = uy/dy;
  float uz = pxyz[2]; int Z = uz/dz;
  if(prlev) {
    printf("%2d ",chx);
    printf("%.3lf %.3lf %.3lf\n",ux,uy,uz);
  }
  if(prlev > 1) {
    printf(" Starting chx=%c (%.3lf,%.3lf,%.3lf)",chx, ux,uy,uz);
    printf(" [%d %d %d]\n",X,Y,Z);
    printf(" dx=(%.3lf,%.3lf,%.3lf) ",dx,dy,dz);
  }

  Ray ray;
  ray.here = -1;
  ray.next = -1;
  ray.dist = INF;
  if((X < 0) || (X >= nx) || (Y < 0) || (Y >= ny) || (Z < 0) || (Z >= nz)) return ray;
  char chh = uv[(Z*ny+Y)*nx+X];
  if(chh != chx) {
    printf("*** %s:%d ",__FILE__,__LINE__);
    printf("  pxyz=(%.3lf,%.3lf,%.3lf)", pxyz[0],pxyz[1],pxyz[2]);
    printf("  dxyz=(%.3lf,%.3lf,%.3lf)",dx,dy,dz);
    printf("  ux=(%.3lf,%.3lf,%.3lf)",ux,uy,uz);
    printf("\n");
    printf("  XYZ=%d %d %d\n",X,Y,Z);
    printf("*** mismatch look chx=%c actual=%c\n",chx,chh);
  }
  assert(chh == chx);
  ray.here = chh;

  float wx = wxyz[0];  float wy = wxyz[1];  float wz = wxyz[2];
  if(prlev > 1) printf(" X=(%d,%d,%d) chx=%c  wx=(%.4lf,%.4lf,%.4lf)\n",X,Y,Z, chx, wx,wy,wz);

  float tMaxX = INF;
  int stepX = 1;
  if(wx > eps) {
    tMaxX = ((X+1)*dx-ux)/wx;
    stepX = 1;
  }
  else if(wx < -eps) {
    tMaxX = (X*dx-ux)/wx;
    stepX = -1;
  }
  float tDeltaX = INF;
  if(fabs(wx) > eps)tDeltaX = stepX*dx/wx;
  
  float tMaxY = INF;
  int stepY = 1;
  if(wy > eps) {
    tMaxY = ((Y+1)*dy-uy)/wy;
    stepY = 1;
  }
  else if(wy < -eps) {
    tMaxY = (Y*dy-uy)/wy;
    stepY = -1;
  }
  float tDeltaY = INF;
  if(fabs(wy) > eps) tDeltaY = stepY*dy/wy;
  float tMaxZ = INF;
  int stepZ = 1;
  if(wz > eps) {
    tMaxZ = ((Z+1)*dz-uz)/wz;
    stepZ = 1;
  }
  else if(wz < -eps) {
    tMaxZ = (Z*dz-uz)/wz;
    stepZ = -1;
  }
  float tDeltaZ = INF;
  if(fabs(wz) > eps) tDeltaZ = stepZ*dz/wz;
  if(prlev > 1) {
   printf(" X=(%d,%d,%d) ",X,Y,Z);
   printf(" chx=%c  U=(%.4lf,%.4lf,%.4lf)",chx, ux,uy,uz);
   printf("\n");
   printf(" tDelta=(%.4lf,%.4lf,%.4lf)",tDeltaX,tDeltaY,tDeltaZ);
   printf(" tMax=(%.4lf,%.4lf,%.4lf)",tMaxX,tMaxY,tMaxZ);
   printf("\n");
  }
  float qx = ux+stepX*tMaxX*wx;
  float qy = uy+stepY*tMaxY*wy;
  float qz = uz+stepZ*tMaxZ*wz;
  if(prlev > 1) printf(" qx=(%.4lf,%.4lf,%.4lf)",qx,qy,qz);
  float ti = tMaxX;
  if(tMaxY < ti) ti = tMaxY;
  if(tMaxZ < ti) ti = tMaxZ;
  if(prlev > 1) {
    printf("  S=(%.4lf,%.4lf,%.4lf)",ux+ti*wx,uy+ti*wy,uz+ti*wz);
    printf("\n");
  }

  /*		*/
  /*  J.Amandatides and A.Woo, A Fast Voxel Traversal Algorithm for Ray Tracing. */
  char chn = chx;
  float t;
  int out = 0;
  while(1) {
    if(tMaxX < tMaxY) {
      if(tMaxX  < tMaxZ) {
        X += stepX;
        tMaxX += tDeltaX;
      } else {
        Z += stepZ;
        tMaxZ += tDeltaZ;
      }
    }
    else {
      if(tMaxY  < tMaxZ) {
        Y += stepY;
        tMaxY += tDeltaY;
      } else {
        Z += stepZ;
        tMaxZ += tDeltaZ;
      }
    }
    if(X < 0)   { out=1; break; }
    if(X >= nx) { out=1; break; }
    if(Y < 0)   { out=1; break; }
    if(Y >= ny) { out=1; break; }
    if(Z < 0)   { out=1; break; }
    if(Z >= nz) { out=1; break; }
    if(prlev > 1) printf(" XYZ=(%d,%d,%d) ",X,Y,Z);
    float qx = ux+tMaxX*wx;
    float qy = uy+tMaxY*wy;
    float qz = uz+tMaxZ*wz;
    if(prlev > 1) printf(" Q=(%.3lf,%.3lf,%.3lf)", qx,qy,qz);
    chn = uv[(Z*ny+Y)*nx+X];
    if(prlev > 1) printf("  chn=%c",chn);
    float tx = tMaxX;
    float ty = tMaxY;
    float tz = tMaxZ;
    t = tx;
    if(ty < t) t = ty;
    if(tz < t) t = tz;
    if(prlev > 1) {
      printf(" t=%.3lf S=(%.4lf,%.4lf,%.4lf)",t, ux+t*wx,uy+t*wy,uz+t*wz);
      printf("\n");
    }
    if(chn != chx) break;
    if(isnan(t)) {
      printf("*** %s:%d t=%.3le\n",__FILE__,__LINE__,t);
      printf("  X=%d Y=%d Z=%d ",X,Y,Z);
      printf("  tMaxX-%.2lf tMaxY=%.2lf tMaxZ=%.2lf\n",
	tMaxX,tMaxY,tMaxZ);
      exit(0);
    }
    if(prlev) {
      printf("%2d ",chx);
      printf("%.3lf %.3lf %.3lf\n",ux+t*wx,uy+t*wy,uz+t*wz);
    }
  }
  if(prlev > 1) printf("  t=%.3lf\n",t);
  if(out) ray.next = -1;
  else ray.next = chn;
  ray.dist = t;
  return ray;
};	// Vox::track
