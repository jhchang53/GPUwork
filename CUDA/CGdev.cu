/*
	CGdev.cu
*/
#include <stdio.h>
#include <assert.h>
#include "CGdev.h"

extern __device__ int dev_ncgs;
extern __device__ int *dev_CGtype;
extern __device__ double *dev_CGparm;

__device__ int trackRPP(double param[], double P[], double trk[]);
__device__ int trackCYL(double param[], double P[], double trk[]);
__device__ int trackTRC(double param[], double P[], double trk[]);

__device__ Ray trackCG(double P[])
{
  int prlev = 0;
  if(prlev) {
    printf(" dev_ncgs=%d\n",dev_ncgs);
    for(int n=0; n < dev_ncgs; n++) {
      printf(" n=%d CGtype=%d CGparm=%.2lf P=(%.3lf,%.3lf,%.3lf)\n",
      n,dev_CGtype[n],dev_CGparm[n*CGBAR],
      P[2],P[3],P[4]);
    }
  }
  /*  find track length of all CG  */
  int *isok = new int[dev_ncgs];
  double *t_min = new double[dev_ncgs];
  double *t_max = new double[dev_ncgs];

  int thiszone = 0;
  double dist = 1.0e+20;
  int zprod = 1;
  double trk[2];

  for(int n=0; n < dev_ncgs; n++) {
    int ok;
    switch(dev_CGtype[n]) {
     case 1:
       ok = trackRPP(dev_CGparm+n*CGBAR,P, trk);
       printf(" RPP ok=%d trk=%.2le %.2le\n",ok,trk[0],trk[1]);
       break;
     case 2:
       ok = trackCYL(dev_CGparm+n*CGBAR,P, trk);
       printf(" CYL ok=%d trk=%.2le %.2le\n",ok,trk[0],trk[1]);
       break;
     case 3:
       ok = trackTRC(dev_CGparm+n*CGBAR,P, trk);
       printf(" TRC ok=%d trk=%.2le %.2le\n",ok,trk[0],trk[1]);
       break;
     default:
       assert(1==0);
    }   // case
    isok[n] = ok;
    if(ok) {
      t_min[n] = trk[0];
      t_max[n] = trk[1];
      /*  find this zone */
      if((trk[0] <= 0.0) && (trk[1] > 0.0)) {
        thiszone += zprod;
      }
      if((t_min[n] > 0.0) && (dist > t_min[n])) dist = t_min[n];
      if((t_max[n] > 0.0) && (dist > t_max[n])) dist = t_max[n];
    }
    zprod = 2*zprod;    // prepare for next CG
  }     // for
  /* check next zone  */
  zprod = 1;
  int nextzone = 0;
  for(int n=0; n < dev_ncgs; n++) {
    if(isok[n]) {
      if((t_min[n] <= dist) && (dist < t_max[n])) {
        nextzone += zprod;
      }
    }
    zprod = 2*zprod;
  }
  printf(" trackCGtrack thiszone=%d nextzone=%d dist=%.5lf\n",thiszone,nextzone,dist);
  Ray ray;
  ray.here = thiszone;
  ray.next = nextzone;
  ray.dist = dist;
  return ray;
};      // trackCG


/*
        CGrpp.cpp
        ref) P.J.Schneider and D.H. Eberly, Geometric Tools for Computer Graphics,
         Elsevier, 2003.
*/

__device__ int Clip(double denom, double numer, double &t0, double &t1)
{
  // Liang{Barsky clipping for a linear component against a face of a box.
  if(denom > 0.0) {
    if(numer > denom*t1) return 0;
    if(numer > denom*t0) t0 = numer/denom;
    return 1;
  }
  else if(denom < 0.0) {
    if(numer > denom*t0) return 0;
    if(numer > denom*t1) t1 = numer/denom;
    return 1;
  }
  else {
    if(numer <= 0.0) return 1;
    else return 0;
  }
  return 0;
};      // CGrpp::Clip


__device__ int trackRPP(double cgparm[], double P8[], double trk[])
{
  double box_C[3],box_e[3];
  for(int i=0; i < 3; i++) {
    box_C[i] = cgparm[2*i];
    box_e[i] = cgparm[2*i+1];
  }
  double INF = 1.0/0.0;
  /*  same as Box::DoLineQuery  */
  double P[3],w[3];
  for(int i=0; i < 3; i++) P[i] = P8[2+i]-box_C[i];
  for(int i=0; i < 3; i++) w[i] = P8[5+i];
  double t0 = -INF;
  double t1 = +INF;
  int numPoints;
  int notCulled =
    Clip(+w[0], -P[0]-box_e[0], t0,t1) &&
    Clip(-w[0], +P[0]-box_e[0], t0,t1) &&
    Clip(+w[1], -P[1]-box_e[1], t0,t1) &&
    Clip(-w[1], +P[1]-box_e[1], t0,t1) &&
    Clip(+w[2], -P[2]-box_e[2], t0,t1) &&
    Clip(-w[2], +P[2]-box_e[2], t0,t1);
  if(notCulled) {
    if(t1 > t0) {
      // The intersection is a segment P + t * w with t in [ t0 , t1 ] .
      numPoints = 2;
      trk[0] = t0;
      trk[1] = t1;
    }
    else {
      // The intersection is a segment P + t * w with t = t0 .
      numPoints = 1;
      trk[0] = t0;
      trk[1] = t0;
    }
  }
  else {
    // The line does not intersect the box. Return invalid parameters .
    numPoints = 0;
    trk[0] = +INF;
    trk[1] = -INF;
  }
  if(numPoints < 2) numPoints = 0;      // fix
  return numPoints;
};	// trackRPP

/*
        ref) P.J.Schneider and D.H. Eberly, Geometric Tools for Computer Graphics,
         Elsevier, 2003.
*/
__device__ int trackCYL(double cgparm[], double P[], double trk[])
{
  int prlev = 0;
  double R = cgparm[0];
  double H0 = cgparm[1];
  double H1 = cgparm[2];
  if(prlev) printf(" R=%.2lf H0=%.2lf H1=%.2lf\n",R,H0,H1);
  /* Linear Components and Cylinders    */
  /* ref. Schneider, chap.11.3.4        */
  double INF = 1.0e+20;
  trk[0] = INF;
  trk[1] = -INF;
  double wx = P[5];  double wy = P[6];  double wz = P[7];
  double a = wx*wx + wy*wy;
  double px = P[2];     double py = P[3];       double pz = P[4];
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
    t[0] = t0;  valid[0] = 1; t[1] = t1;  valid[1] = 1;
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
    double x2 = px+t2*wx;       double y2 = py+t2*wy;
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
};

/*
        P.J.Schneider and D.H.Eerly, Geometric Tools for computer graphics,
        Elsevier (2003)
        chap 11.3.5. Linear components and a cone
*/
__device__ int trackTRC(double cgparm[], double P[], double trk[])
{
  int prlev = 0;
  double Z0 = cgparm[0];
  double Z1 = cgparm[1];
  double Vz = cgparm[2];
  // double H = cgparm[3];
  double gamsq = cgparm[4];

  double INF = 1.0/0.0;
  double eps = 1.0e-20;
  trk[0] = INF;  trk[1] = -INF;
  double wx = P[5];  double wy = P[6];  double wz = P[7];
  if(prlev) printf(" w=(%.4lf,%.4lf,%.4lf)\n",wx,wy,wz);
  double px = P[2];  double py = P[3];  double pz = P[4];
  /*  we know that TRC aligned in z-axis result diagonal M matrix */
  double Mw0 = -gamsq*wx; double Mw1 = -gamsq*wy; double Mw2 = (1-gamsq)*wz;
  double Md0 = -gamsq*px; double Md1 = -gamsq*py; double Md2 = (1-gamsq)*(pz-Vz);
  double c0 = px*Md0+py*Md1+(pz-Vz)*Md2;
  double c1 = wx*Md0+wy*Md1+wz*Md2;
  double c2 = wx*Mw0+wy*Mw1+wz*Mw2;
  if(prlev) printf("c0=%.3le c1=%.3le c2=%.3le",c0,c1,c2);
  double det = c1*c1-c0*c2;
  if(prlev) printf(" det=%.3le\n",det);
  double tmin,tmax;
  if(det > eps) {
    if(fabs(c2) > 1.0e-5) {
      double rtdet = sqrt(det);
      double t0 = (-c1+rtdet)/c2;
      double t0z = pz + t0*wz;
      if(prlev) printf("  t0=%.4lf t0z=%.4lf\n",t0,t0z);
      double t1 = (-c1-rtdet)/c2;
      double t1z = pz + t1*wz;
      if(prlev) printf("  t1=%.4lf t1z=%.4lf\n",t1,t1z);
      /*  take t1 as upper and t0 lower */
      if(t1z < t0z) {
        double texz = t1z;      double tex = t1;
        t1z = t0z;              t1 = t0;
        t0z = texz;             t0 = tex;
      }
      /*  classify  */
      if(Z1 < t0z) {
        // case 1) Vz < Z0 < Z1 < t0z < t1z  : No point
        return 0;
      }
      if(t1z < Vz) {
        // case 10)  t0z < t1z < Vz  < Z0 < Z1 : No point
        return 0;
      }
      if(Z1 < t1z) {    // 2), 4), 7)
        if(Z0 < t0z) {
          // case 2) Vz < Z0 < t0z < Z1 < t1z : t0 and Z1 intersection
          if(prlev) {
            printf(" case 2:\n");
            printf("  Vz=%.3lf Z0=%.3lf t0z=%.3lf Z1=%.3lf t1z=%.3lf\n",
                Vz , Z0 , t0z , Z1 , t1z );
            assert(Vz < Z0);
            assert(Z0 < t0z);
            assert(t0z < Z1);
            assert(Z1 < t1z);
          }
          tmin = t0;
          tmax = (Z1-pz)/wz;
          if(prlev) printf(" == z0=%.3lf zmin=%.3lf(%.2lf) zmax=%.3lf(%.2lf)\n",
                pz, pz+wz*tmin,tmin, pz+wz*tmax,tmax);
        }
        else {
          if(Vz < t0z) {
            // case 4) Vz < t0z < Z0 < Z1 < t1z : Z0 inter and Z1 inter
            if(prlev) {
              printf(" case 4:\n");
              printf("  Vz=%.3lf t0z=%.3lf Z0=%.3lf Z1=%.3lf t1z=%.3lf\n",
                Vz , t0z , Z0 , Z1 , t1z);
              assert(Vz < t0z);
              assert(t0z < Z0);
              assert(Z0 < Z1);
              assert(Z1 < t1z);
            };
            tmin = (Z0-pz)/wz;
            tmax = (Z1-pz)/wz;
          }
          else {
            // case 7) case 7) t0z < Vz < Z0 < Z1 < t1z :  No point
            return 0;
          }
        }
      }
      else {    // 3) 5) 6) 8) 9)
        if(Vz < t0z) {  // 3), 5), 6)
          if(Z0 < t0z) {
            // case 3) Vz < Z0 < t0z < t1z < Z1 : t0 and t1
            if(prlev) {
              printf(" case 3:\n");
              printf(" Vz=%.2lf Z0=%.2lf t0z=%.2lf t1z=%.2lf Z1=%.2lf\n",
                Vz , Z0 , t0z , t1z , Z1);
              assert(Vz < Z0);
              assert(Z0 < t0z);
              assert(t0z < t1z);
              assert(t1z < Z1);
            }
            tmin = t0;
            tmax = t1;
          }
          else {        // 5) 6)
            if(Z0 < t1z) {
              // case 5)  Vz < t0z < Z0 < t1z < Z1 : Z0 inter and t1
              if(prlev) {
                printf(" case 5:\n");
                printf(" Vz=%.2lf t0z=%.2lf Z0=%.2lf t1z=%.2lf Z1=%.2lf\n",
                Vz , t0z , Z0 , t1z , Z1);
              }
              tmin = (Z0-pz)/wz;
              tmax = t1;
            }
            else {
              // case 6)  Vz < t0z < t1z < Z0 < Z1 : No point
              return 0;
            }
          }
        }
        else {  // 8), 9)
          if(Z0 < t1z) {
            // case 8)  t0z < Vz < Z0 < t1z < Z1 : t1 and Z1 inter
            if(prlev) {
              printf(" case 8:");
              printf(" t0z=%.2lf  Vz=%.2lf Z0=%.2lf t1z=%.2lf Z1=%.2lf\n",
                   t0z, Vz, Z0, t1z, Z1);
            }
            tmin = t1;
            tmax = (Z1-pz)/wz;
          }
          else {
            // case 9) t0z < Vz < t1z < Z0 < Z1 : Z0 inter and Z1 inter
            if(prlev) {
              printf(" case 9:\n");
              printf(" t0z=%.2lf Vz=%.2lf t1z=%.2lf Z0=%.2lf Z1=%.2lf\n",
                 t0z, Vz, t1z, Z0, Z1);
            }
            tmin = (Z0-pz)/wz;
            tmax = (Z1-pz)/wz;
          }
        }
      }
      // end of 10 cases
    }   // end of if c2 > 0
    else {
      // parallel to cone surface
      double t = -c0/(2*c1);
      double tz = pz + wz*t;
      if(Z0 < tz) { // case 11) or 12)
        if(Z1 < tz) {
          // case 11: Vz < Z0 < Z1 < tz : no point
          return 0;
        }
        else {
          // case 12: Vz < Z0 < tz < t1 : tz and Z1 inter
          if(prlev) {
            printf(" case 12:\n");
          };

          tmin = t;
          tmax = (Z1-pz)/wz;
        }
      }
      else {    // case 13) or 14)
        if(Vz < tz) {
          // case 13:  Vz < tz < Z0 < Z1 : Z0 inter amd Z1 inter
          if(prlev) {
            printf(" case 13:\n");
          }
          tmin = (Z0-pz)/wz;
          tmax = (Z1-pz)/wz;
        }
        else {
          // case 14:  tz < Vz < Z0 < Z1 : no point
          return 0;
        }
      }
    }   // if c2 = 0
  }     // if det > 0
  else {        // if det = 0
    return 0;
  }
  /*  make tmin < tmax  */
  if(tmin > tmax) {
    double tex = tmin;
    tmin = tmax;
    tmax = tex;
  }
  trk[0] = tmin;
  trk[1] = tmax;
  return 1;
};	// trackTRC


