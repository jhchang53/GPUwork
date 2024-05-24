/*
	SRCdev.cu
*/
#include <stdio.h>
#include <assert.h>
#include <curand_kernel.h>
#include <curand.h>

#include "SRCdev.h"

extern __device__ curandStateMRG32k3a *state;

extern __device__ double dev_Rsrc;
extern __device__ int dev_nsrcg;
extern __device__ double *dev_SRCen;
extern __device__ double *dev_SRCcdf;


__device__ void disk_source(double P[])
{
  int prlev = 0;
  int id = threadIdx.x  + blockIdx.x * blockDim.x;
  curandStateMRG32k3a localState = state[id];

  P[0] = 1.0;	// particle weight
  /*  determine enrgey  using Wheeler algorithm */
  /*  energy in decreasing order */
  double ranen = curand_uniform_double(&localState);
  int ig;
  double weit;
  if(ranen <= dev_SRCcdf[0]) {
    ig = 0;
    weit =  ranen/dev_SRCcdf[0];
  }
  else {
    int ifrom = 0;
    int ito = dev_nsrcg;
    while(ito-ifrom > 1) {
      int iguess = ifrom + (ito-ifrom)/2;
      if(ranen < dev_SRCcdf[iguess]) ito = iguess;
      else ifrom = iguess;
    }
    ig = ito;
    weit = (ranen-dev_SRCcdf[ig-1])/(dev_SRCcdf[ig]-dev_SRCcdf[ig-1]);
    printf(" ig=%d ran=%.5lf cdl=%.5lf cdu=%.5lf weit=%.4lf\n",
      ig,ranen,dev_SRCcdf[ig-1],dev_SRCcdf[ig],weit);
    assert(ig > 0);
  }
  double Ener = weit*dev_SRCen[ig]+(1.0-weit)*dev_SRCen[ig+1];
  printf("ig=%d weit=%.5lf Ener=%.3le\n",ig,weit,Ener);
  P[1] = Ener;
  /*  determine position  */
  double ranx,rany;
  double x,y;
  while(1) {
    ranx = curand_uniform_double(&localState);
    rany = curand_uniform_double(&localState);
    printf("  dev_Rsrc=%.3lf ranx=%.4lf rany=%.4lf\n",dev_Rsrc,ranx,rany);
    x = 2*ranx-1.0;
    y = 2*rany-1.0;
    if(x*x+y*y <= 1.0) break;
  }
  P[2] = dev_Rsrc*x;
  P[3] = dev_Rsrc*y;
  P[4] = 0.0;	// z-position
  printf(" Psrc=(%.4lf,%.4lf)\n",P[2],P[4]);
  
  /*  select angle cosine from equiprobable bin  */
  //  it NEED to be MODIFIED to use ig dependent eqbin
  double eqbin[] = {1.0, 9.42800e-1, 8.80380e-1, 8.15039e-1, 7.47901e-1,
             6.77881e-1, 6.03572e-1, 5.24954e-1, 4.39224e-1, 3.50185e-1,
             2.43306e-1, 1.16176e-1,-3.55876e-2,-1.65392e-1,-2.95196e-1,
            -4.25000e-1,-5.47748e-1,-6.60841e-1,-7.73934e-1,-8.87008e-1,-1.00000e+00};
  int nbin = sizeof(eqbin)/sizeof(double) - 1;
  double rana = curand_uniform_double(&localState);
  int ja = (int) nbin*rana;
  assert(ja < nbin);
  double dxmu = eqbin[ja]-eqbin[ja+1];
  double dprob = 1.0/nbin;
  double xmu = eqbin[ja] + dxmu*(ja*dprob-rana);
  assert(xmu >= -1.0);
  assert(xmu <= 1.0);
  /*  azimuthal angle  using Neumann method, Lux p.21  */
  double rho1,rho2;
  while(1) {
    rho1 = 2*curand_uniform_double(&localState) - 1.0;
    rho2 = 2*curand_uniform_double(&localState) - 1.0;
    if(rho1*rho1+rho2*rho2 <= 1.0) break;
  }
  double sinmu = sqrt(1.0-xmu*xmu);
  double rhosq = rho1*rho1+rho2*rho2;
  double wx = sinmu*(rho1*rho1-rho2*rho2)/rhosq;
  double wy = sinmu*2*rho1*rho2/rhosq;
  if(prlev) {
  double wsq = xmu*xmu + wx*wx + wy*wy;
  printf("  rana=%.4lf ja=%d thup=%.4lf thdn=%.4lf  w=%.4lf %.4lf %.4lf  werr=%.2le\n",
   rana,ja,eqbin[ja],eqbin[ja+1],xmu,wx,wy, 1.0-wsq);
  }
  P[5] = wx;
  P[6] = wy;
  P[7] = xmu;

  state[id] = localState;
};

