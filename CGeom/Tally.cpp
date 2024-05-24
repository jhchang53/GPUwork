/*
	Tally.cpp
	requires atomic operation to accumulate tally
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "Tally.h"

Tally::Tally()
{
  prlev = 0;
  INF = 1.0e+20;
  eps = 1.0e-20;
  nvoxel = 0;
  tal = NULL;
};

Tally::~Tally()
{
  delete [] tal;
};

int Tally::setLattice(const int nxyz[], const double dxyz[], char *uv_i)
{
  nx = nxyz[0];  dx = dxyz[0];
  ny = nxyz[1];  dy = dxyz[1];
  nz = nxyz[2];  dz = dxyz[2];
  printf(" Box: %.3lf,%.3lf,%.3lf\n",nx*dx,ny*dy,nz*dz);
  uv = uv_i;
  /* reserve tally array */
  nvoxel = nz*ny*nx;
  // tal = new float[nvoxel];
  tal = new float[nvoxel];
  for(int ijk=0; ijk < nvoxel; ijk++) tal[ijk] = 0.0;
  return nvoxel;
};

void Tally::clear()
{
  for(int ijk=0; ijk < nvoxel; ijk++) tal[ijk] = 0.0;
};

void Tally::TLE(const double pxyz[], const double wxyz[], const char chx, const double dist)
{
  if(chx <= 32) return;
  if(dist <= 1.0e-10) return;
  /* we are inside UV vox now  */
  double ux = pxyz[0]; int X = ux/dx;
  double uy = pxyz[1]; int Y = uy/dy;
  double uz = pxyz[2]; int Z = uz/dz;
  if((X < 0) || (X >= nx)) {
    printf("*** %s:%d X=%d nx=%d ",__FILE__,__LINE__,X,nx);
    printf("  pxyz[0]=%.3lf ux=%.3lf dx=%.3lf dist=%.3lf\n",pxyz[0],ux,dx, dist);

  }
  assert((X >= 0) && (X < nx));
  assert((Y >= 0) && (Y < ny));
  assert((Z >= 0) && (Z < nz));

  assert(chx == uv[(Z*ny+Y)*nx+X]);
  double wx = wxyz[0];  double wy = wxyz[1];  double wz = wxyz[2];
  if(prlev) {
   printf("  ux=(%.3lf,%.3lf,%.3lf)",ux,uy,uz);
   printf(" XYZ=(%d,%d,%d)",X,Y,Z);
   printf("  w=(%.1lf,%.1lf,%.1lf)\n",wx,wy,wz);
  }
  /*  tally position  */
  int Xt = X;  int Yt = Y; int Zt = Z;
  
  double tMaxX = INF;
  int stepX = 1;
  if(wx > eps) {
    tMaxX = ((X+1)*dx-ux)/wx;
    stepX = 1;
  }
  else if(wx < -eps) {
    tMaxX = (X*dx-ux)/wx;
    stepX = -1;
  }
  double tDeltaX = INF;
  if(fabs(wx) > eps)tDeltaX = stepX*dx/wx;

  double tMaxY = INF;
  int stepY = 1;
  if(wy > eps) {
    tMaxY = ((Y+1)*dy-uy)/wy;
    stepY = 1;
  }
  else if(wy < -eps) {
    tMaxY = (Y*dy-uy)/wy;
    stepY = -1;
  }
  double tDeltaY = INF;
  if(fabs(wy) > eps) tDeltaY = stepY*dy/wy;

  double tMaxZ = INF;
  int stepZ = 1;
  if(wz > eps) {
    tMaxZ = ((Z+1)*dz-uz)/wz;
    stepZ = 1;
  }
  else if(wz < -eps) {
    tMaxZ = (Z*dz-uz)/wz;
    stepZ = -1;
  }
  double tDeltaZ = INF;
  if(fabs(wz) > eps) tDeltaZ = stepZ*dz/wz;

  if(prlev > 1) {
   printf(" X=(%d,%d,%d) ",X,Y,Z);
   printf(" chx=%c  U=(%.4lf,%.4lf,%.4lf)",chx, ux,uy,uz);
   printf("\n");
   printf(" tDelta=(%.4lf,%.4lf,%.4lf)",tDeltaX,tDeltaY,tDeltaZ);
   printf(" tMax=(%.4lf,%.4lf,%.4lf)",tMaxX,tMaxY,tMaxZ);
   printf("\n");
  }
  // double qx = ux+stepX*tMaxX*wx;
  // double qy = uy+stepY*tMaxY*wy;
  // double qz = uz+stepZ*tMaxZ*wz;
  double ti = tMaxX;
  if(tMaxY < ti) ti = tMaxY;
  if(tMaxZ < ti) ti = tMaxZ;
  if(prlev) {
  if(ti >= dist) {
    printf("*** %s:%d ti=%.5le dist=%.5le",__FILE__,__LINE__,ti,dist);
    char chq = uv[(Z*ny+Y)*nx+X];
    printf(" chx=%c  chq=%c",chx,chq);
    printf("  XYZ=(%d,%d,%d)\n",X,Y,Z);
  }
  }
  int out = 0;
  if(ti >= dist) {
    ti = dist;
    out = 1;
  }
  tal[(Zt*ny+Yt)*nx+Xt] += ti;
  if(out) return;
  /*            */
  /*  J.Amandatides and A.Woo, A Fast Voxel Traversal Algorithm for Ray Tracing. */
  char chn = chx;
  double t=0;
  double told = ti;
  // while(t <= dist) {
  while(1) {
    Xt = X;  Yt = Y; Zt = Z;	// tally position
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
    if(X < 0)   break;
    if(X >= nx) break;
    if(Y < 0)   break;
    if(Y >= ny) break;
    if(Z < 0)   break;
    if(Z >= nz) break;
    // qx = ux+tMaxX*wx;
    // qy = uy+tMaxY*wy;
    // qz = uz+tMaxZ*wz;
    chn = uv[(Zt*ny+Yt)*nx+Xt];
    if(chn != chx) {
      printf(" %s:%d chn=%c chx=%c  XYZ=(%d,%d,%d) ",__FILE__,__LINE__,
	chn,chx, Xt,Yt,Zt);
      printf(" dist=%.4lf told=%.4lf",dist,told);
      printf("\n");
    }
    assert(chn == chx);
    double tx = tMaxX;
    double ty = tMaxY;
    double tz = tMaxZ;
    t = tx;
    if(ty < t) t = ty;
    if(tz < t) t = tz;
    if(t >= dist) {
      t = dist;
      out =1;
    }
    assert(t >= told);
    double dtl = t-told;
    told = t;
    if(prlev) printf("  xyz=(%d,%d,%d) dtl=%.3le\n",X,Y,Z, dtl);
    tal[(Zt*ny+Yt)*nx+Xt] += dtl;
    if(prlev > 1) {
      printf(" t=%.3lf S=(%.4lf,%.4lf,%.4lf)",t, ux+t*wx,uy+t*wy,uz+t*wz);
      printf("\n");
    }
    if(out) break;
    // if(chn != chx) break;
    if(prlev > 1) {
      printf("%2d ",chx);
      printf("%.3lf %.3lf %.3lf\n",ux+t*wx,uy+t*wy,uz+t*wz);
    }
  }
  /*  add remaining tally  */
  if(told < dist) {
    double dtl = dist-told;
    assert((Xt >= 0) && (Xt < nx));
    assert((Yt >= 0) && (Yt < ny));
    assert((Zt >= 0) && (Zt < nz));
    tal[(Zt*ny+Yt)*nx+Xt] += dtl;
    told = dist;
  }
  assert(fabs(told-dist) < eps);
};	// Tally::TLE

void Tally::print()
{
  print(tal);
};

void Tally::print(const float txyz[])
{
  for(int i=0; i < nx; i++) {
    int nonz = 0;
    for(int j=0;!nonz && (j < ny); j++) {
      for(int k=0; !nonz && (k < nz); k++) if(txyz[(k*ny+j)*nx+i] > 0) nonz = 1;
    }
    if(nonz) {
      printf("%3d",i);
      for(int k=0; k < nz; k++) {
        double v = 0;
        for(int j=0; j < ny; j++) v += txyz[(k*ny+j)*nx+i];
        if(v > 1.0) printf("*");
        else if(v > 0.1) printf("=");
        else if(v > 0.01) printf("+");
        else if(v > 0.0) printf(".");
        else printf(" ");
      }; // for k
      printf("|\n");
    }
  }
};	// Tally::print
