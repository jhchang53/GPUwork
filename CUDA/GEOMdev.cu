/*
	GEOMdev.cu
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "GEOMdev.h"
#include "CGdev.h"
#include "BOX.h"

extern __device__ int UVnx,UVny,UVnz;
extern __device__ double UVdx,UVdy,UVdz;
extern __device__ char *dev_uv;


__device__ double Rmat[9];
__device__ double Zb,PTx,PTy,PTz;

__device__ void TRAN_set(double target[], double Zb_i, double theta, double phi, double psi)
{
  int prlev = 0;
  Zb = Zb_i;
  PTx = target[0];  PTy = target[1];  PTz = target[2];
  printf(" Zb=%.2lf PT=(%.3lf,%.3lf,%.3lf)\n",Zb,PTx,PTy,PTz);
  /*  rotate z axis by azimutal theta, then polar phi  */
  printf("  polar=%.3lf  azimuth=%.3lf  psi=%.3lf\n",theta,phi,psi);
  double PI = 2.0*acos(0.0);
  double deg_rad = PI/180.0;
  /*  Euler rotation on z axis  */
  double Rz[9];
  for(int ij=0; ij < 9; ij++) Rz[ij] = 0.0;
  double cosphi = cos(deg_rad*phi);
  double sinphi = sin(deg_rad*phi);
  Rz[0] =  cosphi;   Rz[1] = sinphi;
  Rz[3] = -sinphi;   Rz[4] = cosphi;
  Rz[8] = 1.0;
  /*  Euler rotation on y axis  */
  double Ry[9];
  for(int ij=0; ij < 9; ij++) Ry[ij] = 0.0;
  double costheta = cos(deg_rad*theta);
  double sintheta = sin(deg_rad*theta);
  Ry[0] =  costheta;  Ry[2] = -sintheta;
  Ry[4] = 1.0;
  Ry[6] = sintheta;  Ry[8] = costheta;

  /*  find Rodrigous vector for rotation  */
  double ki[] = {0.0,0.0,1.0};
  double kiz[3];
  for(int j=0; j < 3; j++) {
    double sum = 0.0;
    for(int i=0; i < 3; i++) sum += Rz[i*3+j]*ki[i];
    kiz[j] = sum;
  }
  printf("kiz: (%.3lf,%.3lf,%.3lf)\n",kiz[0],kiz[1],kiz[2]);

  double kizy[3];
  for(int j=0; j < 3; j++) {
    double sum = 0.0;
    for(int i=0; i < 3; i++) sum += Ry[i*3+j]*kiz[i];
    kizy[j] = sum;
  }
  printf("kizy: (%.3lf,%.3lf,%.3lf)\n",kizy[0],kizy[1],kizy[2]);
  double kx =  cosphi*sintheta;
  double ky = -sinphi*sintheta;
  double kz =  costheta;
  printf(" k=(%.3lf,%.3lf,%.3lf)\n",kx,ky,kz);
  double K[9];
  K[0] =   0;  K[1] = -kz;  K[2] =  ky;
  K[3] =  kz;  K[4] = 0;    K[5] = -kx;
  K[6] = -ky;  K[7] =  kx;  K[8] = 0.0;
  for(int i=0; i < 3; i++) {
    printf("K %d:",i);
    for(int j=0; j < 3; j++) printf(" %.2le",K[i*3+j]);
    printf("\n");
  }

  double K2[9];
  for(int i=0; i < 3; i++) {
    for(int k=0; k < 3; k++) {
      double sum = 0.0;
      for(int j=0; j < 3; j++) {
        sum += K[i*3+j]*K[j*3+k];
      }
      K2[i*3+k] = sum;
    }
  }
  for(int i=0; i < 3; i++) {
    printf("Ksq %d:",i);
    for(int j=0; j < 3; j++) printf(" %.2le",K2[i*3+j]);
    printf("\n");
  }
  /*  Rodrigues' rotation matrix  */
  double Rr[9];
  /*  R = I + sinpsi*K + (1-cospsi)*K^2  */
  for(int ij=0; ij < 9; ij++) Rr[ij] = 0.0;
  for(int i=0; i < 3; i++) Rr[i*3+i] = 1.0;
  double sinpsi = sin(deg_rad*psi);
  double cospsi = cos(deg_rad*psi);
  for(int ij=0; ij < 9; ij++) Rr[ij] += sinpsi*K[ij]+(1.0-cospsi)*K2[ij];
  for(int i=0; i < 3; i++) {
    printf("Rr %d:",i);
    for(int j=0; j < 3; j++) printf(" %.2le",Rr[i*3+j]);
    printf("\n");
  }
  /*  compose full rotation matrix  */
  double Ryz[9];
  for(int i=0; i < 3; i++) {
    for(int j=0; j < 3; j++) {
      double sum = 0.0;
      for(int k=0; k < 3; k++) {
        sum += Ry[i*3+k]*Rz[j*3+k];
      }
      Ryz[i*3+j] = sum;
    }
  }
  for(int i=0; i < 3; i++) {
    for(int j=0; j < 3; j++) {
      double sum = 0.0;
      for(int k=0; k < 3; k++) {
        sum += Rr[i*3+k]*Ryz[j*3+k];
      }
      Rmat[i*3+j] = sum;
    }
  }
  for(int i=0; i < 3; i++) {
    printf("R  %d:",i);
    for(int j=0; j < 3; j++) printf(" %.2le",Rmat[i*3+j]);
    printf("\n");
  }
  if(prlev > 1) {
    /*  check orthogonality  */
    double O[9];
    for(int i=0; i < 3; i++) {
      for(int k=0; k < 3; k++) {
        double sum = 0.0;
        for(int j=0; j < 3; j++) {
          sum +=Rmat[i*3+j]*Rmat[k*3+j];
        }
        O[i*3+k] = sum;
      }
    }
    for(int i=0; i < 3; i++) {
      printf("O %d:",i);
      for(int j=0; j < 3; j++) printf(" %.4lf",O[i*3+j]);
      printf("\n");
    }
    assert(1==0);
  }
};	// TRAN_set

__device__ void TRAN_S2U(double Q[], double P[])
{
  /*  convert point in source geometry into point of UV voxel */
  int prlev = 0;
  P[0] = Q[0];	// weight
  P[1] = Q[1];	// energy
  if(prlev) printf(" Q=(%.3lf,%.3lf,%.3lf)\n",Q[2],Q[3],Q[4]);

  /*  vector from target (to be rotated  */
  double pb[3];
  pb[0] = - Q[2];	// offset 2 for particle weight and energy
  pb[1] = - Q[3];
  pb[2] = Zb - Q[4];
  if(prlev) printf(" Pbar=(%.3lf,%.3lf,%.3lf)\n",pb[0],pb[1],pb[2]);
  /*  rotated vector  */
  double Rp[3];
  for(int i=0; i < 3; i++) {
    double sum = 0;
    for(int j=0; j < 3; j++) {
      sum += Rmat[i*3+j]*pb[j];
    }
    Rp[i] = sum;
  }
  /*  vector from origin  */
  double px = Rp[0] + PTx;
  double py = Rp[1] + PTy;
  double pz = Rp[2] + PTz;
  if(prlev) printf(" P=(%.3lf,%.3lf,%.3lf)\n",px,py,pz);
  P[2] = px;
  P[3] = py;
  P[4] = pz;
  /*  convert vector in source geometry into vector of UV voxel */
  /*  vector from target (to be rotated  */
  double pv[3];
  pv[0] = -Q[5];
  pv[1] = -Q[6];
  pv[2] = -Q[7];
  if(prlev) printf(" Pvec=(%.3lf,%.3lf,%.3lf)\n",pv[0],pv[1],pv[2]);
  /*  rotated vector  */
  for(int i=0; i < 3; i++) {
    double sum = 0;
    for(int j=0; j < 3; j++) {
      sum += Rmat[i*3+j]*pv[j];
    }
    P[5+i] = sum;
  }
};	// TRAN_S2U

__device__ void TRAN_U2S(double P[], double Q[])
{
  /* conver point in UV voxel geometry into point in source geom  */
  /*  vector from point of rotation  */
  double pt[3];
  pt[0] = P[2] - PTx;
  pt[1] = P[3] - PTy;
  pt[2] = P[4] - PTz;
  /*  inverse rotated */
  double IR[3];
  for(int i=0; i < 3; i++) {
    double sum = 0;
    for(int j=0; j < 3; j++) {
      sum += Rmat[j*3+i]*pt[j];    // for orthogonal matrix inv R is transpose
    }
    IR[i] = sum;
  }

  double qx = -IR[0];  double qy = -IR[1];  double qz = Zb-IR[2];
  Q[2] = qx;  Q[3] = qy;  Q[4] = qz;
  /*  direction vector from point of rotation  */
  double pv[3];
  pv[0] = P[5];  pv[1] = P[6];  pv[2] = P[7];
  for(int i=0; i < 3; i++) {
    double sum = 0;
    for(int j=0; j < 3; j++) {
      sum += Rmat[j*3+i]*pv[j]; 	// for orthogonal matrix inv R is transpose
    }
    Q[5+i] = -sum;
  }
};      // TRAN_U2S


__device__ char getChar(double P[])
{
  /*  retrieve char at pxyz, if outside returns 0  */
  int X = P[2]/UVdx;
  if((X < 0) || (X >= UVnx)) return 0;
  int Y = P[3]/UVdy;
  if((Y < 0) || (Y >= UVny)) return 0;
  int Z = P[4]/UVdz;
  if((Z < 0) || (Z >= UVnz)) return 0;
  return dev_uv[(Z*UVny+Y)*UVnx+X];
};	// getChar

__device__ Ray trackUV(double P[], char chx)
{
  int prlev = 0;
  /* orgin coordinate of lattice is (0,0,0)  */
  double ux = P[2]; int X = ux/UVdx;
  double uy = P[3]; int Y = uy/UVdy;
  double uz = P[4]; int Z = uz/UVdz;
  if(prlev) {
    printf("%2d ",chx);
    printf("%.3lf %.3lf %.3lf\n",ux,uy,uz);
  }
  if(prlev > 1) {
    printf(" Starting chx=%c (%.3lf,%.3lf,%.3lf)",chx, ux,uy,uz);
    printf(" [%d %d %d]\n",X,Y,Z);
    printf(" dx=(%.3lf,%.3lf,%.3lf) ",UVdx,UVdy,UVdz);
  }

  double INF = 1.0/0.0;
  Ray ray;
  ray.here = -1;
  ray.next = -1;
  ray.dist = INF;
  if((X < 0) || (X >= UVnx) || (Y < 0) || (Y >= UVny) || (Z < 0) || (Z >= UVnz)) return ray;
  char chh = dev_uv[(Z*UVny+Y)*UVnx+X];
  if(chh != chx) {
    printf("*** %s:%d ",__FILE__,__LINE__);
    printf("  pxyz=(%.3lf,%.3lf,%.3lf)", P[2],P[3],P[4]);
    printf("  dxyz=(%.3lf,%.3lf,%.3lf)",UVdx,UVdy,UVdz);
    printf("  ux=(%.3lf,%.3lf,%.3lf)",ux,uy,uz);
    printf("\n");
    printf("  XYZ=%d %d %d\n",X,Y,Z);
    printf("*** mismatch look chx=%c actual=%c\n",chx,chh);
  }
  assert(chh == chx);
  ray.here = chh;

  double wx = P[5];  double wy = P[6];  double wz = P[7];
  if(prlev > 1) printf(" X=(%d,%d,%d) chx=%c  wx=(%.4lf,%.4lf,%.4lf)\n",X,Y,Z, chx, wx,wy,wz);

  double tMaxX = INF;
  double eps = 1.0e-20;
  int stepX = 1;
  if(wx > eps) {
    tMaxX = ((X+1)*UVdx-ux)/wx;
    stepX = 1;
  }
  else if(wx < -eps) {
    tMaxX = (X*UVdx-ux)/wx;
    stepX = -1;
  }
  double tDeltaX = INF;
  if(fabs(wx) > eps)tDeltaX = stepX*UVdx/wx;

  double tMaxY = INF;
  int stepY = 1;
  if(wy > eps) {
    tMaxY = ((Y+1)*UVdy-uy)/wy;
    stepY = 1;
  }
  else if(wy < -eps) {
    tMaxY = (Y*UVdy-uy)/wy;
    stepY = -1;
  }
  double tDeltaY = INF;
  if(fabs(wy) > eps) tDeltaY = stepY*UVdy/wy;
  double tMaxZ = INF;
  int stepZ = 1;
  if(wz > eps) {
    tMaxZ = ((Z+1)*UVdz-uz)/wz;
    stepZ = 1;
  }
  else if(wz < -eps) {
    tMaxZ = (Z*UVdz-uz)/wz;
    stepZ = -1;
  }
  double tDeltaZ = INF;
  if(fabs(wz) > eps) tDeltaZ = stepZ*UVdz/wz;
  if(prlev > 1) {
   printf(" X=(%d,%d,%d) ",X,Y,Z);
   printf(" chx=%c  U=(%.4lf,%.4lf,%.4lf)",chx, ux,uy,uz);
   printf("\n");
   printf(" tDelta=(%.4lf,%.4lf,%.4lf)",tDeltaX,tDeltaY,tDeltaZ);
   printf(" tMax=(%.4lf,%.4lf,%.4lf)",tMaxX,tMaxY,tMaxZ);
   printf("\n");
  }
  double qx = ux+stepX*tMaxX*wx;
  double qy = uy+stepY*tMaxY*wy;
  double qz = uz+stepZ*tMaxZ*wz;
  if(prlev > 1) printf(" qx=(%.4lf,%.4lf,%.4lf)",qx,qy,qz);
  double ti = tMaxX;
  if(tMaxY < ti) ti = tMaxY;
  if(tMaxZ < ti) ti = tMaxZ;
  if(prlev > 1) {
    printf("  S=(%.4lf,%.4lf,%.4lf)",ux+ti*wx,uy+ti*wy,uz+ti*wz);
    printf("\n");
  }

  /*            */
  /*  J.Amandatides and A.Woo, A Fast Voxel Traversal Algorithm for Ray Tracing. */
  char chn = chx;
  double t;
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
    if(X >= UVnx) { out=1; break; }
    if(Y < 0)   { out=1; break; }
    if(Y >= UVny) { out=1; break; }
    if(Z < 0)   { out=1; break; }
    if(Z >= UVnz) { out=1; break; }
    if(prlev > 1) printf(" XYZ=(%d,%d,%d) ",X,Y,Z);
    double qx = ux+tMaxX*wx;
    double qy = uy+tMaxY*wy;
    double qz = uz+tMaxZ*wz;
    if(prlev > 1) printf(" Q=(%.3lf,%.3lf,%.3lf)", qx,qy,qz);
    chn = dev_uv[(Z*UVny+Y)*UVnx+X];
    if(prlev > 1) printf("  chn=%c",chn);
    double tx = tMaxX;
    double ty = tMaxY;
    double tz = tMaxZ;
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
      assert(1==0);
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
};      // trackUV

__device__ double CGP[8];

__device__ Ray howfar(double P[])
{
  int prlev = 1;
  double INF = 1.0/0.0;
  Ray rc;
  rc.dist = INF;
  /*  p, w is UV  */
  char chx = getChar(P);
  int inUV  = 0;
  if(chx >= 32) inUV = 1;
  if(inUV) {
    if(prlev) printf("== %s:%d inside of UV\n",__FILE__,__LINE__);
    if(chx < 32) chx = 32;
    // rc = vox->track(p,w,chx);
    rc = trackUV(P,chx);
    if(prlev) {
      printf("  Vox here=%c next=%c dist=%.2le\n",chx,rc.next, rc.dist);
    }
    chx = rc.next;
    if(chx < 0) {
      inUV = 0;
      return rc;
    }
    // printf(" at line %d  chx=%c dist=%.2le\n",__LINE__,chx,rc.dist);
    if(chx != ' ')  return rc;
    // printf("  passing line %d\n",__LINE__);
    /*  next cell is void  */
    if(prlev) printf("== %s:%d next cell is void\n",__FILE__,__LINE__);
    /*  check all CG if there is intersection */
    
    // trc->U2S(p,CGp);
    // trc->U2Svec(w,CGw);
    // Ray cgray = howfarCG(CGp,CGw);
    TRAN_U2S(P,CGP);
    Ray cgray = trackCG(CGP);
    if(prlev) printf("CG  ray=%d %d %.3lf\n",cgray.here,cgray.next,cgray.dist);
    // int isout = vox->checkXYZ();
    if(cgray.here == 0) return rc;
    if(prlev) printf("== %s:%d next cell is void\n",__FILE__,__LINE__);
    rc.next = -1;
    if(prlev) printf("**  escape\n");
    inUV = 0;
    // printf("  passing line %d\n",__LINE__);
    rc.dist = 0.0;
    return rc;
  }
  else {        // not in a UV box
    /* it can be CG or UV  */
    // printf("  passing line %d\n",__LINE__);
    // trc->U2S(p,CGp);
    // trc->U2Svec(w,CGw);
    // Ray cgray = howfarCG(CGp,CGw);
    TRAN_U2S(P,CGP);
    Ray cgray = trackCG(CGP);
    if((cgray.here == 0) && (cgray.dist > 1.0e+10)) {
      // printf("  passing line %d\n",__LINE__);
      // check UV box */
      double t[2];
      int meetBox = BOX_DoRayQuery(P, t);
      if(prlev) printf("Box  meet=%d t0=%.3le t1=%.3lf\n",meetBox,t[0],t[1]);
      if(!meetBox) {
        rc.here = -1;  rc.next = -1;
        rc.dist = 0.0;
        if(prlev) printf("**  escape\n");
        return rc;
      }
      else {
        rc.here = 0;
        rc.next = 0;
        rc.dist = t[0];
        if(rc.dist < 0.0) rc.dist = 0;
        inUV = 1;
        if(prlev) printf("#  entering Box at %.2le\n",t[0]);
        return rc;
      }
    }
    else {
      if(cgray.dist < 1.0e-10) cgray.dist = 1.0e-10;
      return cgray;
    }
  }
  printf("== %s:%d\n",__FILE__,__LINE__);
  assert(1==0);
  return rc;
};	// howfar
