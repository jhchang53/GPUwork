/*
   	CGhost.cu
*/
#include "CGhost.h"

/*  type 1 */
int CGrpp_host(double *cgparm, const double a[], const double b[])
{
  int cgtype = 1;
  double box_C[3],box_e[3];
  for(int i=0; i < 3; i++) {
    box_C[i] = (b[i]+a[i])/2;
    box_e[i] = (b[i]-a[i])/2;
  }
  for(int i=0; i < 3; i++) {
    cgparm[2*i] = box_C[i];
    cgparm[2*i+1] = box_e[i];
  }
  return cgtype;
};	//  CGrpp_host

/*  type 2 */
int CGcyl_host(double *cgparm,  const double R,
  const double H0, const double H1)
{
  int cgtype = 2;
  cgparm[0] = R;
  cgparm[1] = H0;
  cgparm[2] = H1;
  return cgtype;
};

/*  type 3  */
int CGtrc_host(double *cgparm, const double R0, const double R1,
  const double Z0, const double Z1)
{
  int cgtype = 3;
  double Vz = (R1*Z0-R0*Z1)/(R1-R0);
  double R0sq = R0*R0;
  double R1sq = R1*R1;
  double H = Z1-Z0;
  double gamsq = H*H/((R1-R0)*(R1-R0)+H*H);
  cgparm[0] = Z0;
  cgparm[1] = Z1;
  cgparm[2] = Vz;
  cgparm[3] = H;
  cgparm[4] = gamsq;

  return cgtype;
};	// CGtrc_host
