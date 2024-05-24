/*
	Voxel.cpp
	ref. J.Amandatides and A.Woo, A Fast Voxel Traversal Algorithm
	for Ray Tracing.

*/
#include <stdio.h>
#include <math.h>
#include "Voxel.h"

Voxel::Voxel()
{
};

Voxel::~Voxel()
{
};

void Voxel::setLat2D(double dx_i, double x0_i, int nx_i,
  double dy_i, double y0_i, int ny_i, char*  uv_i)
{
  dx = dx_i;
  x0 = x0_i;
  nx = nx_i;
  xb = x0+nx*dx;
  dy = dx;
  y0 = y0_i;
  ny = ny_i;
  yb = y0+ny*dy;
  uv = uv_i;
};	// Voxel::setLat2D

int Voxel::traverse2D(double px, double py, double wx, double wy)
{
  /*  check whether (px,py) is inside of 2D lattice */
  if(px < x0) return 0;
  if(px >= xb) return 0;
  if(py < y0) return 0;
  if(py >= yb) return 0;
  double ux = px-x0;
  int X = ux/dx;
  double uy = py-y0;
  int Y = uy/dy;
  char chx = uv[Y*nx+X];
  printf(" dx=(%.3lf,%.3lf) ",dx,dy);
  printf(" X=%d Y=%d chx=%c  wx=(%.4lf,%.4lf)",X,Y, chx, wx,wy);
  double eps = 1.0e-20;
  double tMaxX = 1.0e+20;
  int stepX = 1;
  if(wx > eps) {
    tMaxX = ((X+1)*dx-ux)/wx;
    stepX = 1;
  }
  else if(wx < eps) {
    tMaxX = (X*dx-ux)/wx;
    stepX = -1;
  }
  double tDeltaX = stepX*dx/wx;
  double tMaxY = 1.0e+20;
  int stepY = 1;
  if(wy > eps) {
    tMaxY = ((Y+1)*dy-uy)/wy;
    stepY = 1;
  }
  else if(wy < eps) {
    tMaxY = (Y*dy-uy)/wy;
    stepY = -1;
  }
  double tDeltaY = stepY*dy/wy;
  printf(" X=%d Y=%d chx=%c  U=(%.4lf,%.4lf)",X,Y, chx, ux,uy);
  printf("\n");
  printf(" tDelta=(%.4lf,%.4lf)",tDeltaX,tDeltaY);

  printf(" tMax=(%.4lf,%.4lf)",tMaxX,tMaxY);
  double qx = ux+stepX*tMaxX*wx;
  double qy = uy+stepY*tMaxY*wy;
  printf(" qx=(%.4lf,%.4lf)",qx,qy);
  double ti = tMaxX;
  if(tMaxY < ti) ti = tMaxY;
  printf("  S=(%.4lf,%.4lf)",ux+ti*wx,uy+ti*wy);
  printf("\n");
  /*  fast traverse  */
  // while((X >= 0) && (X < nx-1) && (Y >= 0) && (Y < ny-1)) {
  char chi = chx;
  while(1) {
    if(tMaxX < tMaxY) {
      tMaxX += tDeltaX;
      X += stepX;
      if(X < 0) break;
      if(X >= nx) break;
    }
    else {
      tMaxY += tDeltaY;
      Y += stepY;
      if(Y < 0) break;
      if(Y >= ny) break;
    }
    printf(" XY=(%d,%d) tMax=(%.4lf,%.4lf)",X,Y, tMaxX,tMaxY);
    double qx = ux+tMaxX*wx;
    double qy = uy+tMaxY*wy;
    printf(" Q=(%.3lf,%.3lf)", qx,qy);
    chi = uv[Y*nx+X];
    printf("  ch=%c",chi);
    double tx = tMaxX;
    double ty = tMaxY;
    printf("  tx=(%.4lf,%.4lf)",tx,ty);
    double t = tx;
    if(ty < t) t = ty;
    printf("  S=(%.4lf,%.4lf)",ux+t*wx,uy+t*wy);
    printf("\n");
  }
  return 1;
};	// Voxel::traverse2D

void Voxel::setLat3D(double dxyz[], double xyz0[], int nxyz[], char *uv_i)
{
  dx = dxyz[0];   dy = dxyz[1];   dz = dxyz[2];
  x0 = xyz0[0];   y0 = xyz0[1];   z0 = xyz0[2];
  nx = nxyz[0];   ny = nxyz[1];   nz = nxyz[2];
  xb = x0+nx*dx;  yb = y0+ny*dy;  zb = z0+nz*dz;
  uv = uv_i;
};	// Voxel::setLat3D

int Voxel::traverse3D(double pxyz[], double wxyz[])
{
  /*  check whether (px,py) is inside of 2D lattice */
  double px = pxyz[0];
  if(px < x0) return 0;
  if(px >= xb) return 0;
  double py = pxyz[1];
  if(py < y0) return 0;
  if(py >= yb) return 0;
  double pz = pxyz[2];
  if(pz < z0) return 0;
  if(pz >= zb) return 0;
  double ux = px-x0;
  int X = ux/dx;
  double uy = py-y0;
  int Y = uy/dy;
  double uz = pz-z0;
  int Z = uz/dz;

  char chx = uv[(Z*ny+Y)*nx+X];
  printf(" dx=(%.3lf,%.3lf,%.3lf) ",dx,dy,dz);
  double wx = wxyz[0];  double wy = wxyz[1];  double wz = wxyz[2];
  printf(" X=(%d,%d,%d) chx=%c  wx=(%.4lf,%.4lf,%.4lf)\n",X,Y,Z, chx, wx,wy,wz);
  double eps = 1.0e-20;
  double tMaxX = 1.0e+20;
  int stepX = 1;
  if(wx > eps) {
    tMaxX = ((X+1)*dx-ux)/wx;
    stepX = 1;
  }
  else if(wx < eps) {
    tMaxX = (X*dx-ux)/wx;
    stepX = -1;
  }
  double tDeltaX = stepX*dx/wx;
  double tMaxY = 1.0e+20;
  int stepY = 1;
  if(wy > eps) {
    tMaxY = ((Y+1)*dy-uy)/wy;
    stepY = 1;
  }
  else if(wy < eps) {
    tMaxY = (Y*dy-uy)/wy;
    stepY = -1;
  }
  double tDeltaY = stepY*dy/wy;
  double tMaxZ = 1.0e+20;
  int stepZ = 1;
  if(wz > eps) {
    tMaxZ = ((Z+1)*dz-uz)/wz;
    stepZ = 1;
  }
  else if(wz < eps) {
    tMaxZ = (Z*dz-uz)/wz;
    stepZ = -1;
  }
  double tDeltaZ = stepZ*dz/wz;
  printf(" X=(%d,%d,%d) ",X,Y,Z);
  printf(" chx=%c  U=(%.4lf,%.4lf,%.4lf)",chx, ux,uy,uz);
  printf("\n");
  printf(" tDelta=(%.4lf,%.4lf,%.4lf)",tDeltaX,tDeltaY,tDeltaZ);

  printf(" tMax=(%.4lf,%.4lf,%.4lf)",tMaxX,tMaxY,tMaxZ);
  double qx = ux+stepX*tMaxX*wx;
  double qy = uy+stepY*tMaxY*wy;
  double qz = uz+stepZ*tMaxZ*wz;
  printf(" qx=(%.4lf,%.4lf,%.4lf)",qx,qy,qz);
  double ti = tMaxX;
  if(tMaxY < ti) ti = tMaxY;
  if(tMaxZ < ti) ti = tMaxZ;
  printf("  S=(%.4lf,%.4lf,%.4lf)",ux+ti*wx,uy+ti*wy,uz+ti*wz);
  printf("\n");
  /*  fast traverse  */
  char chi = chx;
  while(1) {
    if(tMaxX < tMaxY) {
      if(tMaxX  < tMaxZ) {
        X += stepX;
        if(X < 0) break;
        if(X >= nx) break;
        tMaxX += tDeltaX;
      } else {
        Z += stepZ;
        if(Z < 0) break;
        if(Z >= nz) break;
        tMaxZ += tDeltaZ;
      }
    }
    else {
      if(tMaxY  < tMaxZ) {
        Y += stepY;
        if(Y < 0) break;
        if(Y >= ny) break;
        tMaxY += tDeltaY;
      } else {
        Z += stepZ;
        if(Z < 0) break;
        if(Z >= nz) break;
        tMaxZ += tDeltaZ;
      }
    }
    printf(" XYZ=(%d,%d,%d) tMax=(%.4lf,%.4lf,%.4lf)",X,Y,Z, tMaxX,tMaxY,tMaxZ);
    double qx = ux+tMaxX*wx;
    double qy = uy+tMaxY*wy;
    double qz = uz+tMaxZ*wz;
    printf(" Q=(%.3lf,%.3lf,%.3lf)", qx,qy,qz);
    chi = uv[(Z*ny+Y)*nx+X];
    if(chi != ' ') printf("  ch=%c",chi);
    double tx = tMaxX;
    double ty = tMaxY;
    double tz = tMaxZ;
    if(chi != ' ') printf("  tx=(%.4lf,%.4lf,%.4lf)",tx,ty,tz);
    double t = tx;
    if(ty < t) t = ty;
    if(tz < t) t = tz;
    if(chi != ' ') printf("  S=(%.4lf,%.4lf,%.4lf)",ux+t*wx,uy+t*wy,uz+t*wz);
    printf("\n");
  }
  return 1;
};	//  Voxel::traverse3D

