/*
	CGcyl.cpp
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


void CGcyl::CYL(const double R_i, const double H0_i, const double H1_i)
{
  R = R_i;
  assert(R > eps);
  H0 = H0_i;   H1 = H1_i;
  assert(H1 > H0);
};	//  CGcyl::CYL

int CGcyl::track(const double p[], const double w[], double trk[])
{
  /* Linear Components and Cylinders	*/
  /* ref. Schneider, chap.11.3.4  	*/
  trk[0] = INF;
  trk[1] = -INF;
  double wx = w[0];	double wy = w[1];	double wz = w[2];
  double a = wx*wx + wy*wy;
  double px = p[0];	double py = p[1];	double pz = p[2];
  double b = 2*(wx*px+wy*py);
  double c = px*px+py*py-R*R;
  double det = b*b-4*a*c;
  int valid[4];
  double t[4];
  double tmin = INF;
  double tmax = -INF;
  if(det >= 0.0) {
    double rtdet = sqrt(det);
    double t0 = (-b+rtdet)/(2*a);
    double t1 = (-b-rtdet)/(2*a);
    if(prlev) printf(" t=%.3lf %.3lf\n",t0,t1);
    t[0] = t0;	valid[0] = 1; t[1] = t1;  valid[1] = 1;
    double q0x= px+t0*wx;  double q0y = py+t0*wy;  double q0z = pz+t0*wz;
    if((q0z < H0) || (q0z > H1)) valid[0] = 0;
    double q1x= px+t1*wx;  double q1y = py+t1*wy;  double q1z = pz+t1*wz;
    if((q1z < H0) || (q1z > H1)) valid[1] = 0;
    if(prlev) {
      printf("t0=%.3lf t1=%.1lf ",t0,t1);
      printf(" Q0=(%.4lf,%.4lf,%.4lf)[%d]",q0x,q0y,q0z, valid[0]);
      printf(" Q1=(%.4lf,%.4lf,%.4lf)[%d]",q1x,q1y,q1z, valid[1]);
      printf("\n");
    }
    //  check end caps
    double Rsq = R*R;
    double t2 = (H0-pz)/wz;  valid[2] = 1;
    t[2] = t2;
    double x2 = px+t2*wx;	double y2 = py+t2*wy;
    if(x2*x2+y2*y2 > Rsq) valid[2] = 0;
    double q2x= px+t2*wx;  double q2y = py+t2*wy;  double q2z = pz+t2*wz;
    if(prlev) {
      printf(" t2=%.3lf ",t2);
      printf(" Q2=(%.4lf,%.4lf,%.4lf)[%d]\n",q2x,q2y,q2z, valid[2]);
    }
    double t3 = (H1-pz)/wz;  valid[3] = 1;
    t[3] = t3;
    double x3 = px+t3*wx;       double y3 = py+t3*wy;
    if(x3*x3+y3*y3 > Rsq) valid[3] = 0;
    double q3x= px+t3*wx;  double q3y = py+t3*wy;  double q3z = pz+t3*wz;
    if(prlev) {
    printf(" t3=%.3lf ",t3);
    printf(" Q3=(%.4lf,%.4lf,%.4lf)[%d]\n",q3x,q3y,q3z, valid[3]);
    }
    //  select min and max
    tmin = INF;
    tmax = -INF;
    for(int n=0; n < 4; n++) {
      if(valid[n]) {
        double tplus = t[n];
        double qx= px+tplus*wx;  double qy = py+tplus*wy;  double qz = pz+tplus*wz;
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
    double t0 = -b/(2*a);
    double q0x= px+t0*wx;  double q0y = py+t0*wy;  double q0z = pz+t0*wz;
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
