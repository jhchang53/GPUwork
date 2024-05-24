/*
	CGcyl.cpp	= float version
	ref) P.J.Schneider and D.H. Eberly, Geometric Tools for Computer Graphics,
	 Elsevier, 2003. 
*/
#include <stdio.h>
#include <limits>
#include <math.h>
#include <assert.h>
#include "CGcyl.h"


CGcyl::CGcyl()
{
};

CGcyl::~CGcyl()
{
};


void CGcyl::CYL(const float R_i, const float H0_i, const float H1_i)
{
  R = R_i;
  assert(R > eps);
  H0 = H0_i;   H1 = H1_i;
  assert(H1 > H0);
};	//  CGcyl::CYL

int CGcyl::track(const float p[], const float w[], float trk[])
{
  /* Linear Components and Cylinders	*/
  /* ref. Schneider, chap.11.3.4  	*/
  trk[0] = INF;
  trk[1] = -INF;
  float wx = w[0];	float wy = w[1];	float wz = w[2];
  float a = wx*wx + wy*wy;
  float px = p[0];	float py = p[1];	float pz = p[2];
  float b = 2*(wx*px+wy*py);
  float c = px*px+py*py-R*R;
  float det = b*b-4*a*c;
  int valid[4];
  float t[4];
  float tmin = INF;
  float tmax = -INF;
  if(det >= 0.0) {
    float rtdet = sqrt(det);
    float t0 = (-b+rtdet)/(2*a);
    float t1 = (-b-rtdet)/(2*a);
    if(prlev) printf(" t=%.3lf %.3lf\n",t0,t1);
    t[0] = t0;	valid[0] = 1; t[1] = t1;  valid[1] = 1;
    float q0x= px+t0*wx;  float q0y = py+t0*wy;  float q0z = pz+t0*wz;
    if((q0z < H0) || (q0z > H1)) valid[0] = 0;
    float q1x= px+t1*wx;  float q1y = py+t1*wy;  float q1z = pz+t1*wz;
    if((q1z < H0) || (q1z > H1)) valid[1] = 0;
    if(prlev) {
      printf("t0=%.3lf t1=%.1lf ",t0,t1);
      printf(" Q0=(%.4lf,%.4lf,%.4lf)[%d]",q0x,q0y,q0z, valid[0]);
      printf(" Q1=(%.4lf,%.4lf,%.4lf)[%d]",q1x,q1y,q1z, valid[1]);
      printf("\n");
    }
    //  check end caps
    float Rsq = R*R;
    float t2 = (H0-pz)/wz;  valid[2] = 1;
    t[2] = t2;
    float x2 = px+t2*wx;	float y2 = py+t2*wy;
    if(x2*x2+y2*y2 > Rsq) valid[2] = 0;
    float q2x= px+t2*wx;  float q2y = py+t2*wy;  float q2z = pz+t2*wz;
    if(prlev) {
      printf(" t2=%.3lf ",t2);
      printf(" Q2=(%.4lf,%.4lf,%.4lf)[%d]\n",q2x,q2y,q2z, valid[2]);
    }
    float t3 = (H1-pz)/wz;  valid[3] = 1;
    t[3] = t3;
    float x3 = px+t3*wx;       float y3 = py+t3*wy;
    if(x3*x3+y3*y3 > Rsq) valid[3] = 0;
    float q3x= px+t3*wx;  float q3y = py+t3*wy;  float q3z = pz+t3*wz;
    if(prlev) {
    printf(" t3=%.3lf ",t3);
    printf(" Q3=(%.4lf,%.4lf,%.4lf)[%d]\n",q3x,q3y,q3z, valid[3]);
    }
    //  select min and max
    tmin = INF;
    tmax = -INF;
    for(int n=0; n < 4; n++) {
      if(valid[n]) {
        float tplus = t[n];
        float qx= px+tplus*wx;  float qy = py+tplus*wy;  float qz = pz+tplus*wz;
        if(prlev) {
        printf(" tplus=%.3lf ",tplus);
        printf(" Q=(%.4lf,%.4lf,%.4lf)\n",qx,qy,qz);
        }
        if(t[n] < tmin) tmin = t[n];
        if(t[n] > tmax) tmax = t[n];
      }
    }
    if(prlev) printf("  tmin=%.3le tmax=%.3le\n",tmin,tmax);
    if(tmin > tmax) return 0;
  }
  else if(det == 0.0) {
    // ray is tangent to side, no need to check caps
    float t0 = -b/(2*a);
    float q0x= px+t0*wx;  float q0y = py+t0*wy;  float q0z = pz+t0*wz;
    if(q0z > H1) return 0;
    if(q0z < H0) return 0;
    if(prlev) {
    printf("tt=%.3lf ",t0);
    printf(" Qt=(%.4lf,%.4lf,%.4lf)\n",q0x,q0y,q0z);
    }
  }
  else {
   return 0;
  }
  trk[0] = tmin;
  trk[1] = tmax;
  return (tmax > tmin);
};	// CGcyl::CYL
