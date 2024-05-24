/*
	track.cu

*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <curand.h>
#include "SRCdev.h"
#include "CGdev.h"
#include "track.h"
#include "GEOMdev.h"
#include "BOX.h"

__device__ curandStateMRG32k3a *state;

__device__ double dev_Rsrc;
__device__ int dev_nsrcg;
__device__ double *dev_SRCen;
__device__ double *dev_SRCcdf;

__device__ int dev_ncgs;
__device__ int *dev_CGtype;
__device__ double *dev_CGparm;

__device__ int UVnx,UVny,UVnz;
__device__ double UVdx,UVdy,UVdz;
__device__ char *dev_uv;

__device__ float *dev_tally;

/*  initialize GPU variables  */
__global__ void setup_kernel(curandStateMRG32k3a *d_state,
  int ncgs, int *CGtype, double *CGparm,
  int *nxyz, double *dxyz, char *uv, float *tally,
  double Rsrc, int nsrcg, double *SRCen, double *SRCcdf,
  double target[], double Zb, double theta, double phi, double psi)
{
  int prlev = 0;
  /*  initialize random number generator  */
  int id = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  unsigned long long seed = id+293;
  unsigned long long sequence = 1;
  unsigned long long offset = 71;
  state = d_state;
  curand_init(seed,sequence,offset,&state[id]);
  /*  initialize Source  */
  dev_Rsrc = Rsrc;
  dev_nsrcg = nsrcg;
  dev_SRCen = SRCen;
  for(int n=0; n <= nsrcg; n++) dev_SRCen[n] = SRCen[n];
  dev_SRCcdf = SRCcdf;
  for(int n=0; n < nsrcg; n++) dev_SRCcdf[n] = SRCcdf[n];
  if(prlev) {
  printf("SRCen:");
  for(int n=0; n < nsrcg; n++) printf(" %.2le",dev_SRCen[n]);
  printf("\n");
  }
  /*  initialize CG  */
  printf(" Ncgs=%d\n",ncgs);
  dev_ncgs = ncgs;
  printf(" d_ncgs=%d\n",dev_ncgs);
  dev_CGtype = CGtype;
  assert(dev_CGtype != NULL);
  dev_CGparm = CGparm;
  assert(dev_CGparm != NULL);
  for(int n=0; n < ncgs; n++) {
    dev_CGtype[n] = CGtype[n];
    for(int i=0; i < CGBAR; i++) dev_CGparm[n*CGBAR+i] = CGparm[n*CGBAR+i];
  }
  if(prlev) {
    for(int n=0; n < ncgs; n++) {
      printf(" n=%d CGtype=%d parm=%.2lf\n",n,dev_CGtype[n],dev_CGparm[n*CGBAR]);
    }
  }
  /*  initialize UV  */
  UVnx = nxyz[0];   UVny = nxyz[1];  UVnz = nxyz[2];
  UVdx = dxyz[0];   UVdy = dxyz[1];  UVdz = dxyz[2];
  dev_uv = uv;
  dev_tally = tally;
  printf("  nxyz=%d %d %d\n",UVnx,UVny,UVnz);
  printf("  dxyz=%.3lf %.3lf %.3lf\n",UVdx,UVdy,UVdz);
  int tsize = UVnx*UVny*UVnz*sizeof(float);
  /*  set Box for UV region  */
  double bmin[3],bmax[3];
  for(int i=0; i < 3; i++) {
    bmin[i] = 0.0;
    bmax[i] = nxyz[i]*dxyz[i];
  }
  BOX_set(bmin,bmax);

  /*  set coordinate transformation  */
  TRAN_set(target,Zb,theta,phi,psi);
  printf("-- setup -- size=%d\n",tsize);
};

__global__ void generate_kernel()
{
  int prlev = 0;
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  // unsigned int x;
  double x;
  /* Copy state to local memory for efficiency */
  // curandStateMRG32k3a localState = state[id];
  // localState = state[id];
  /* Generate pseudo-random unsigned ints */
  // for(int i = 0; i < 1; i++) {
    // x = curand(&localState);
    // printf(" id=%d x=%d\n",id,x);
    // x = curand_uniform_double(&localState);
    // printf(" id=%d x=%.5lf\n",id,x);
  // }
  double Psrc[8];
  disk_source(Psrc);
  printf(" Psrc=%.3lf %.2le (%.3lf,%.3lf,%.3lf) (%.3lf,%.3lf,%.3lf)\n",
	Psrc[0],Psrc[1],Psrc[2],Psrc[3],Psrc[4],Psrc[5],Psrc[6],Psrc[7]);
  Ray ray = trackCG(Psrc);
  printf("Ray: %d %d %.2le\n",ray.here,ray.next,ray.dist);
  if(prlev) {
    for(int n=0; n < dev_ncgs; n++) {
      printf(" n=%d CGtype=%d \n",n,dev_CGtype[n]); // d_CGtype[n]);
    }
  }
  
  float talval = 1.0;
  int jx = threadIdx.x;  int jy = blockIdx.x; int jz = 1;
  atomicAdd(dev_tally+(jz*UVny+jy)*UVnx+jx, talval);
};
