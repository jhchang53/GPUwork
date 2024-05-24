/*
	TRcoord.cpp - float version
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TRcoord.h"

TRcoord::TRcoord()
{
  prlev = 0;
  PI = 2*acos(0.0);
  deg_rad = PI/180.0;
};

TRcoord::~TRcoord()
{
};

void TRcoord::setPrlev(int prl)
{
  prlev = prl;
};

void TRcoord::set(const float target[], const float Zb_i, const float theta, const float phi,
  const float psi)
{
  Zb = Zb_i;
  PTx = target[0];  PTy = target[1];  PTz = target[2];
  printf(" Zb=%.2lf PT=(%.3lf,%.3lf,%.3lf)\n",Zb,PTx,PTy,PTz);
  /*  rotate z axis by azimutal theta, then polar phi  */
  printf("  polar=%.3lf  azimuth=%.3lf  psi=%.3lf\n",theta,phi,psi);
  /*  Euler rotation on z axis  */
  float Rz[9];
  for(int ij=0; ij < 9; ij++) Rz[ij] = 0.0;
  float cosphi = cos(deg_rad*phi);
  float sinphi = sin(deg_rad*phi);
  Rz[0] =  cosphi;   Rz[1] = sinphi;
  Rz[3] = -sinphi;   Rz[4] = cosphi;
  Rz[8] = 1.0;
  /*  Euler rotation on y axis  */
  float Ry[9];
  for(int ij=0; ij < 9; ij++) Ry[ij] = 0.0;
  float costheta = cos(deg_rad*theta);
  float sintheta = sin(deg_rad*theta);
  Ry[0] =  costheta;  Ry[2] = -sintheta;
  Ry[4] = 1.0;
  Ry[6] = sintheta;  Ry[8] = costheta;

  /*  find Rodrigous vector for rotation  */
  float ki[] = {0.0,0.0,1.0};
  float kiz[3];
  for(int j=0; j < 3; j++) {
    float sum = 0.0;
    for(int i=0; i < 3; i++) sum += Rz[i*3+j]*ki[i];
    kiz[j] = sum;
  }
  printf("kiz: (%.3lf,%.3lf,%.3lf)\n",kiz[0],kiz[1],kiz[2]);

  float kizy[3];
  for(int j=0; j < 3; j++) {
    float sum = 0.0;
    for(int i=0; i < 3; i++) sum += Ry[i*3+j]*kiz[i];
    kizy[j] = sum;
  }
  printf("kizy: (%.3lf,%.3lf,%.3lf)\n",kizy[0],kizy[1],kizy[2]);


  /*  axis of rotation (theta,phi) */
  /*  theta is downward  */
  // float costheta = cos(deg_rad*theta);	// downward
  // float sintheta = sin(deg_rad*theta);
 //  float cosphi = cos(deg_rad*phi);
  //  float sinphi = sin(deg_rad*phi);
  float kx =  cosphi*sintheta;
  float ky = -sinphi*sintheta;
  float kz =  costheta;
  printf(" k=(%.3lf,%.3lf,%.3lf)\n",kx,ky,kz);
  float K[9];
  K[0] =   0;  K[1] = -kz;  K[2] =  ky;
  K[3] =  kz;  K[4] = 0;    K[5] = -kx;
  K[6] = -ky;  K[7] =  kx;  K[8] = 0.0;
  for(int i=0; i < 3; i++) {
    printf("K %d:",i);
    for(int j=0; j < 3; j++) printf(" %.2le",K[i*3+j]);
    printf("\n");
  }

  float K2[9];
  for(int i=0; i < 3; i++) {
    for(int k=0; k < 3; k++) {
      float sum = 0.0;
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
  float Rr[9];
  /*  R = I + sinpsi*K + (1-cospsi)*K^2  */
  for(int ij=0; ij < 9; ij++) Rr[ij] = 0.0;
  for(int i=0; i < 3; i++) Rr[i*3+i] = 1.0;
  float sinpsi = sin(deg_rad*psi);
  float cospsi = cos(deg_rad*psi);
  for(int ij=0; ij < 9; ij++) Rr[ij] += sinpsi*K[ij]+(1.0-cospsi)*K2[ij];
  for(int i=0; i < 3; i++) {
    printf("Rr %d:",i);
    for(int j=0; j < 3; j++) printf(" %.2le",Rr[i*3+j]);
    printf("\n");
  }
  /*  compose full rotation matrix  */
  float Ryz[9];
  for(int i=0; i < 3; i++) {
    for(int j=0; j < 3; j++) {
      float sum = 0.0;
      for(int k=0; k < 3; k++) {
        sum += Ry[i*3+k]*Rz[j*3+k];
      }
      Ryz[i*3+j] = sum;
    }
  }
  for(int i=0; i < 3; i++) {
    for(int j=0; j < 3; j++) {
      float sum = 0.0;
      for(int k=0; k < 3; k++) {
        sum += Rr[i*3+k]*Ryz[j*3+k];
      }
      R[i*3+j] = sum;
    }
  }
  for(int i=0; i < 3; i++) {
    printf("R  %d:",i);
    for(int j=0; j < 3; j++) printf(" %.2le",R[i*3+j]);
    printf("\n");
  }
  if(prlev > 1) {
    /*  check orthogonality  */
    float O[9];
    for(int i=0; i < 3; i++) {
      for(int k=0; k < 3; k++) {
        float sum = 0.0;
        for(int j=0; j < 3; j++) {
          sum += R[i*3+j]*R[k*3+j];
        }
        O[i*3+k] = sum;
      }
    }
    for(int i=0; i < 3; i++) {
      printf("O %d:",i);
      for(int j=0; j < 3; j++) printf(" %.4lf",O[i*3+j]);
      printf("\n");
    }
    exit(0);
  }
  /*  translation */
  Tx = target[0] + Zb*cosphi*sintheta;
  Ty = target[1] - Zb*sinphi*sintheta;
  Tz = target[2] + Zb*costheta;
  printf(" Tx=%.3lf %.3lf %.3lf\n",Tx,Ty,Tz);
};

void TRcoord::S2U(const float q[], float p[])
{
  /*  convert point in source geometry into point of UV voxel */
  if(prlev) printf(" Q=(%.3lf,%.3lf,%.3lf)\n",q[0],q[1],q[2]);
  /*  vector from target (to be rotated  */
  float pb[3];
  pb[0] = - q[0];
  pb[1] = - q[1];
  pb[2] = Zb - q[2];
  if(prlev) printf(" Pbar=(%.3lf,%.3lf,%.3lf)\n",pb[0],pb[1],pb[2]);
  /*  rotated vector  */
  float Rp[3];
  for(int i=0; i < 3; i++) {
    float sum = 0;
    for(int j=0; j < 3; j++) {
      sum += R[i*3+j]*pb[j];
    }
    Rp[i] = sum;
  }
  /*  vector from origin  */
  float px = Rp[0] + PTx;
  float py = Rp[1] + PTy;
  float pz = Rp[2] + PTz;
  if(prlev) printf(" P=(%.3lf,%.3lf,%.3lf)\n",px,py,pz);
  p[0] = px;
  p[1] = py;
  p[2] = pz;
};	// TRcoord::S2U

void TRcoord::S2Uvec(const float qvec[], float pvec[])
{
  /*  convert vector in source geometry into vector of UV voxel */
  if(prlev) printf(" Q=(%.3lf,%.3lf,%.3lf)\n",qvec[0],qvec[1],qvec[2]);
  /*  vector from target (to be rotated  */
  float pb[3];
  pb[0] = - qvec[0];
  pb[1] = - qvec[1];
  pb[2] = - qvec[2];
  if(prlev) printf(" Pbar=(%.3lf,%.3lf,%.3lf)\n",pb[0],pb[1],pb[2]);
  /*  rotated vector  */
  for(int i=0; i < 3; i++) {
    float sum = 0;
    for(int j=0; j < 3; j++) {
      sum += R[i*3+j]*pb[j];
    }
    pvec[i] = sum;
  }
  if(prlev) printf(" P=(%.3lf,%.3lf,%.3lf)\n",pvec[0],pvec[1],pvec[2]);

};      // TRcoord::S2Uvec

void TRcoord::U2S(const float p[], float q[])
{
  /* conver point in UV voxel geometry into point in source geom  */
  /*  vector from point of rotation  */
  float pt[3];
  pt[0] = p[0] - PTx;
  pt[1] = p[1] - PTy;
  pt[2] = p[2] - PTz;
  /*  inverse rotated */
  float IR[3];
  for(int i=0; i < 3; i++) {
    float sum = 0;
    for(int j=0; j < 3; j++) {
      sum += R[j*3+i]*pt[j];	// for orthogonal matrix inv R is transpose
    }
    IR[i] = sum;
  }

  float qx = -IR[0];  float qy = -IR[1];  float qz = Zb-IR[2];
  if(prlev) printf(" q=(%.3lf,%.3lf,%.3lf)\n",qx,qy,qz);
  q[0] = qx;  q[1] = qy;  q[2] = qz;
};	// TRcoord::U2S

void TRcoord::U2Svec(const float pvec[], float qvec[])
{
  /* conver point in UV voxel geometry into point in source geom  */
  /*  vector from point of rotation  */
  /*  inverse rotated */
  float IR[3];
  for(int i=0; i < 3; i++) {
    float sum = 0;
    for(int j=0; j < 3; j++) {
      sum += R[j*3+i]*pvec[j];    // for orthogonal matrix inv R is transpose
    }
    qvec[i] = -sum;
  }

  if(prlev) printf(" q=(%.3lf,%.3lf,%.3lf)\n",qvec[0],qvec[1],qvec[2]);
};      // TRcoord::U2Svec

