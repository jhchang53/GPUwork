/*
	CGtrc.cpp
	P.J.Schneider and D.H.Eerly, Geometric Tools for computer graphics, 
	Elsevier (2003)
	chap 11.3.5. Linear components and a cone
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "CGtrc.h"

CGtrc::CGtrc()
{
};

CGtrc::~CGtrc()
{
};

void CGtrc::TRC(const double R0_i, const double R1_i, const double Z0_i, const double Z1_i)
{
  R0 = R0_i;   R1 = R1_i;
  Z0 = Z0_i;   Z1 = Z1_i;
  if(prlev) printf("R0=%.4lf R1=%.4lf Z0=%.4lf Z1=%.4lf\n",R0,R1,Z0,Z1);
  Vz = (R1*Z0-R0*Z1)/(R1-R0);
  R0sq = R0*R0;
  R1sq = R1*R1;
  H = Z1-Z0;
  gamsq = H*H/((R1-R0)*(R1-R0)+H*H);
  if(prlev) printf("Vz=%.4lf gamsq=%.4le\n",Vz,gamsq);
};	// CGtrc::TRC

int CGtrc::track(const double p[], const double w[], double trk[])
{
  trk[0] = INF;  trk[1] = -INF;
  double wx = w[0];   double wy = w[1];  double wz = w[2];
  if(prlev) printf(" w=(%.4lf,%.4lf,%.4lf)\n",wx,wy,wz);
  double px = p[0];   double py = p[1];  double pz = p[2];
  /*  we know that TRC aligned in z-axis result diagonal M matrix */
  double Mw0 = -gamsq*wx; double Mw1 = -gamsq*wy; double Mw2 = (1-gamsq)*wz;
  double Md0 = -gamsq*px; double Md1 = -gamsq*py; double Md2 = (1-gamsq)*(pz-Vz);
  double c0 = px*Md0+py*Md1+(pz-Vz)*Md2;
  double c1 = wx*Md0+wy*Md1+wz*Md2;
  double c2 = wx*Mw0+wy*Mw1+wz*Mw2;
  if(prlev) printf("c0=%.3le c1=%.3le c2=%.3le",c0,c1,c2);
  double det = c1*c1-c0*c2;
  if(prlev) printf(" det=%.3le\n",det);
  int kase;
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
        double texz = t1z;	double tex = t1;
        t1z = t0z;		t1 = t0;
        t0z = texz;		t0 = tex;
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
      if(Z1 < t1z) {	// 2), 4), 7)
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
      else {	// 3) 5) 6) 8) 9)
        if(Vz < t0z) {	// 3), 5), 6)
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
          else {	// 5) 6)
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
        else {	// 8), 9)
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
    }	// end of if c2 > 0
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
      else {	// case 13) or 14)
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
    }	// if c2 = 0
  }	// if det > 0
  else {	// if det = 0
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
};	// CGtrc::TRC

