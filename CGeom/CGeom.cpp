/*
	CGeom.cpp
	ref) P.J.Schneider and D.H. Eberly, Geometric Tools for Computer Graphics,
	 Elsevier, 2003. 
*/
#include <stdio.h>
#include <limits>
#include <math.h>
#include "CGeom.h"

using namespace std;

CGeom::CGeom()
{
  prlev = 0;
  epsilon = 1.0e-10;
};

CGeom::~CGeom()
{
};

double CGeom::trackMin()
{
  return tmin;
};

double CGeom::trackMax()
{
  return tmax;
};

int CGeom::RPP(double xmin[], double xmax[], double p[], double w[])
{
  /*  Linear Component and Axis-Aligned Bounding Box	*/
  /*  Schneider, Chap.11.12.2	*/
  double tFar  =  numeric_limits<double>::infinity();
  double tNear =  -tNear;
  for(int n=0; n < 3; n++) {	// for each pair of plnaes
    if(fabs(w[n]) < epsilon) {
      // ray parallel to planes
      if((p[n] < xmin[n]) || (p[n] < xmax[n])) return 0;
    }
    double t0 = (xmin[n]-p[n])/w[n];
    double t1 = (xmax[n]-p[n])/w[n];
    printf("t0=%.3lf t1=%.3lf\n",t0,t1);
    // check ordering
    if(t0 > t1) {	// swap
      double tmp = t1;
      t0 = t1;
      t1 = tmp;
    }
    // compare
    if(t0 > tNear) tNear= t0;
    if(t1 < tFar)  tFar = t1;
    if(tNear > tFar) return 0;
    if(tFar < 0) return 0;
    printf(" tNear=%.3lf tFar=%.3lf\n",tNear,tFar);
  }
  double x0 = p[0] + tNear*w[0];
  double y0 = p[1] + tNear*w[1];
  double z0 = p[2] + tNear*w[2];
  printf("P0=(%.3lf,%.3lf,%.3lf) ",x0,y0,z0);
  double x1 = p[0] + tFar*w[0];
  double y1 = p[1] + tFar*w[1];
  double z1 = p[2] + tFar*w[2];
  printf(" P1=(%.3lf,%.3lf,%.3lf)\n",x1,y1,z1);

  // Box intersect
  double tIntersect;
  if(tNear > 0) tIntersect = tNear;
  else tIntersect = tFar;
  return 1;
};	// CGeom::RPP

int CGeom::CYL(double R, double H0, double H1, double p[], double w[])
{
  /* Linear Components and Cylinders	*/
  /* ref. Schneider, chap.11.3.4  	*/
  double wx = w[0];	double wy = w[1];	double wz = w[2];
  double a = wx*wx + wy*wy;
  double px = p[0];	double py = p[1];	double pz = p[2];
  double b = 2*(wx*px+wy*py);
  double c = px*px+py*py-R*R;
  double det = b*b-4*a*c;
  int valid[4];
  double t[4];
  if(det > 0.0) {
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
    tmin = +1.0e+10;
    tmax = -1.0e+10;
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
    printf("*** %s:%d\n",__FILE__,__LINE__);
    exit(0);
  }
  return 1;
};	// CGeom::RPP

/* 2D lattice with same dx,dy */
void CGeom::setLat2D(double dx_i, double x0_i, int nx_i, double y0_i, int ny_i,
char *uv_i)
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
};	// CGeom::setLat2D

int CGeom::locate2D(double px, double py, double wx, double wy)
{
  /*  check whether (px,py) is inside of 2D lattice */
  if(px < x0) return 0;
  if(px >= xb) return 0;
  if(py < y0) return 0;
  if(py >= yb) return 0;
  ux = px-x0;
  ix = ux/dx;
  uy = py-y0;
  jy = uy/dy;
  char ch = uv[jy*nx+ix];
  printf(" ix=%d jy=%d ch=%c  (%.4lf,%.4lf)\n",ix,jy, ch, px,py);
  int face = -1;
  double sNear;
  if(wx > epsilon) {	// it is not face 2
    sNear = ((ix+1)*dx-ux)/wx;
    face = 0;
    if(wy > epsilon) {
      double s = ((jy+1)*dy-uy)/wy;
      if(s < sNear) {
        sNear = s;
        face = 1;
      }
     }
     else if(wy < -epsilon) {
       double s = (jy*dy-uy)/wy;
       if(s < sNear) {
        sNear = s;
        face = 3;
      }
    }
  }
  else if(wx < -epsilon) {
    sNear = (ix*dx-ux)/wx;
    face = 2;
    if(wy > epsilon) {
      double s = ((jy+1)*dy-uy)/wy;
      if(s < sNear) {
        sNear = s;
        face = 1;
      }
    }
    else if(wy < -epsilon) {
      double s = (jy*dy-uy)/wy;
      if(s < sNear) {
        sNear = s;
        face = 3;
      }
    }
  }
  else {
    if(wy > epsilon) {
      sNear = ((jy+1)*dy-uy)/wy;
      face = 2;
    }
    else if(wy < -epsilon) {
      sNear = (jy*dx-uy)/wy;
      face = 3;
    }
  }
  double fx = ux+sNear*wx;
  double fy = uy+sNear*wy;
  printf("  face=%d sNear=%.3le  (%.4lf,%.4lf)\n",face,sNear,fx,fy);
  for(int jj=0; (jj < 3) && next2D(ch,face,wx,wy,ix,jy); jj++);
  return 1;
};	// CGeom::locate2D

int CGeom::next2D(int ch, int face, double wx, double wy, int ixi, int jyi)
{
  printf("  next2D face=%d  ixi=%d jyi=%d\n",face, ixi,jyi);
  if(face == 0) {	// from left side
    if(ixi+1 >= nx) return 0;
    if(uv[jyi*nx+ixi+1] != ch) {
      printf(" next ch=%c at (%d,%d)\n",uv[jyi*nx+ixi+1], ixi+1,jyi);
      return 0;
    }
    double s = ((ixi+2)*dx-ux)/wx;
    double ycut = uy+s*wy;
    int nextface = 0;
    double sNear = s;
    printf("sx=%.4lf\n",s);
    if(ycut > (jyi+1)*dy) {
      nextface = 1;	jy++;
      sNear = ((jyi+1)*dy-uy)/wy;
    }
    else if(ycut < jyi*dy) {
      nextface = 3;	jy--;
      sNear = (jyi*dy-uy)/wy;
    }
    else {
      ix++;
    }
    double fx = ux+sNear*wx;
    double fy = uy+sNear*wy;
    printf(" nextface=%d sNear=%.3le  (%.4lf,%.4lf)\n",nextface,sNear,fx,fy);
    face = nextface;
  }
  else if(face == 1) {	// from top side of (ixi,jyi)
    if(jyi+1 >= ny) return 0;
    if(uv[(jyi+1)*nx+ixi] != ch) {
      printf(" next ch=%c at (%d,%d)\n",uv[(jyi+1)*nx+ixi], ixi,jyi+1);
      return 0;
    }
    double s = ((jyi+2)*dy-uy)/wy;
    double xcut = ux+wx*s;
    int nextface = 1;
    double sNear = s;
    printf("s=%.4lf\n",s);
    if(xcut >= (ixi+1)*dx) {
      nextface = 0;  ix++;
      sNear = ((ixi+1)*dx-ux)/wx;
    }
    else if(xcut < ixi*dx) {
      nextface = 2;  ix--;
      sNear = (ixi*dx-ux)/wx;
    }
    else {
      jy++;
    }
    double fx = ux+sNear*wx;
    double fy = uy+sNear*wy;
    printf(" nextface=%d sNear=%.3le  (%.4lf,%.4lf)\n",nextface,sNear,fx,fy);
    face = nextface;
  }
  else if(face == 2) {
    if((ixi-1) < 0) return 0;
    if(uv[jyi*nx+ixi-1] != ch) {
      printf(" next ch=%c at (%d,%d)\n",uv[jyi*nx+ixi-1], ixi-1,jyi);
      return 0;
    }
    double s = ((ixi-1)*dx-ux)/wx;
    double ycut = uy+s*wy;
    int nextface = 0;
    double sNear = s;
    printf("sx=%.4lf\n",s);
    if(ycut > (jyi+1)*dy) {
      nextface = 1;     jy++;
      sNear = ((jyi+1)*dy-uy)/wy;
    }
    else if(ycut < jyi*dy) {
      nextface = 3;     jy--;
      sNear = (jyi*dy-uy)/wy;
    }
    else {
      ix--;
    }
    double fx = ux+sNear*wx;
    double fy = uy+sNear*wy;
    printf(" nextface=%d sNear=%.3le  (%.4lf,%.4lf)\n",nextface,sNear,fx,fy);
    face = nextface;
  }
  else {
    printf("** face=%d\n",face);
    exit(0);
  }
  return 1;
};	// CGeom::next2D
