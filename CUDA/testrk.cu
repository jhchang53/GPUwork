/*
	test driver

*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <chrono>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <curand.h>
#include "CGhost.h"
#include "SRChost.h"
#include "track.h"

using namespace std::chrono;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


int *cgtype;
double *cgparm;


/*  define CG geometry  */
int CGinit_host()
{
  int ncgs = 3;
  cgtype = (int *) malloc(ncgs*sizeof(int));
  cgparm = (double *) malloc(ncgs*CGBAR*sizeof(double));
  double xmin[] = {-17.55,-17.55,-0.1};
  double xmax[] = { 17.55, 17.55,10.56};
  cgtype[0] = CGrpp_host(cgparm, xmin,xmax);
  cgtype[1] = CGtrc_host(cgparm+CGBAR, 8.0,17.55,1.0,9.55);
  cgtype[2] = CGcyl_host(cgparm+2*CGBAR, 8.0, 0.0,1.0);
  return ncgs;
};	//  CGinit_host
  
int SRCinit_host();

int main()
{
  system_clock::time_point wtime_begin = system_clock::now();
  const unsigned int threadsPerBlock = THREADS_PER_BLOCK;
  const unsigned int blockCount = BLOCK_COUNT;
  const unsigned int totalThreads = threadsPerBlock * blockCount;
  /*  initialize random number */
  cudaError_t rc;
  curandStateMRG32k3a *devMRGStates;
  rc = cudaMalloc((void **)&devMRGStates, totalThreads *
	sizeof(curandStateMRG32k3a));
  assert(rc == cudaSuccess);	// check GPU memory allocation
  /*  initialize source  */
  int nsrcg = SRCinit_host();

  double Rsrc = 1.0;
  double *d_SRCen;
  rc = cudaMalloc(&d_SRCen,(1+nsrcg)*sizeof(double));
  assert(rc == cudaSuccess);

  cudaMemcpy(d_SRCen,SRCen,(1+nsrcg)*sizeof(double), cudaMemcpyHostToDevice);
  double *d_SRCcdf;
  rc = cudaMalloc(&d_SRCcdf,nsrcg*sizeof(double));
  cudaMemcpy(d_SRCcdf,SRCcdf,nsrcg*sizeof(double), cudaMemcpyHostToDevice);

  
  /*  initialize CG  */
  int ncgs = CGinit_host();
  for(int n=0; n < ncgs; n++) {
    printf(" type=%d cgparm=%.3le %.3le\n",cgtype[n],cgparm[n*CGBAR],
      cgparm[n*CGBAR+1]);
  }
  int *d_CGtype;
  rc = cudaMalloc(&d_CGtype, ncgs*sizeof(int));
  cudaMemcpy(d_CGtype,cgtype,ncgs*sizeof(int), cudaMemcpyHostToDevice);
  double *d_CGparm;
  rc = cudaMalloc(&d_CGparm, ncgs*CGBAR*sizeof(double));
  cudaMemcpy(d_CGparm,cgparm,ncgs*CGBAR*sizeof(double), cudaMemcpyHostToDevice);

  /*  initialize UV  */
  int nxyz[] = {128,128,134};
  double dxyz[] = {0.2,0.2,0.2};
  const char *uvpath = "/home/jhchang/Dicom3/test3/0831150844478.uv";
  FILE *F = fopen(uvpath,"r");
  if(F == NULL) {
    printf("*** %s not exist.\n",uvpath);
    exit(0);
  }
  int nx = nxyz[0]; int ny = nxyz[1];   int nz = nxyz[2];
  int nvoxel = nx*ny*nz;
  char *uv = new char[nvoxel];
  fgets(uv,nvoxel,F);
  fclose(F);
  int *d_nxyz;
  rc = cudaMalloc(&d_nxyz, 3*sizeof(int));
  cudaMemcpy(d_nxyz,nxyz, 3*sizeof(int), cudaMemcpyHostToDevice);
  double *d_dxyz;
  rc = cudaMalloc(&d_dxyz, 3*sizeof(double));
  cudaMemcpy(d_dxyz,dxyz, 3*sizeof(double), cudaMemcpyHostToDevice);
  char *d_uv;
  rc = cudaMalloc(&d_uv, nvoxel);
  cudaMemcpy(d_uv,uv, nvoxel, cudaMemcpyHostToDevice);

  /*  set coord. translation  */
  double Zb = 20.0; // 20.0;
  double theta = 85.0; // 45.0;
  double phi = 0.0;
  double psi = 0.0;
  double target[] = {12.0,12.0,25.0};
  double *d_target;
  rc = cudaMalloc(&d_target, 3*sizeof(double));
  cudaMemcpy(d_target,target, 3*sizeof(double), cudaMemcpyHostToDevice);


  /*  initialize tally  */
  float *d_tally;
  rc = cudaMalloc(&d_tally, nvoxel*sizeof(float));
  cudaMemset(d_tally, 0.0,nvoxel);

  /*  initialize kernel  */
  setup_kernel<<<BLOCK_COUNT,THREADS_PER_BLOCK>>>(devMRGStates,
    ncgs,d_CGtype,d_CGparm,  d_nxyz,d_dxyz,d_uv,  d_tally,
    Rsrc,nsrcg,d_SRCen,d_SRCcdf,
    d_target,Zb,theta,phi,psi);
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
  cudaDeviceSynchronize();

  for(int iter=0; iter < 2; iter++) {
    generate_kernel<<<BLOCK_COUNT,THREADS_PER_BLOCK>>>();
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
    cudaDeviceSynchronize();
    printf(" iter=%d done\n",iter);
  }

  /* collect tally  */
  float *tal = (float *) malloc(nvoxel*sizeof(float));
  cudaMemcpy(tal,d_tally,nvoxel*sizeof(float), cudaMemcpyDeviceToHost);
  int jz = 1;
  for(int jy=0; jy < BLOCK_COUNT; jy++) {
    printf("tal%d:",jy);
    for(int jx=0; jx < THREADS_PER_BLOCK; jx++) printf(" %.2lf",tal[(jz*ny+jy)*nx+jx]);
    printf("\n");
  }
  system_clock::time_point wtime_end = system_clock::now();
  microseconds elapsed_usec = duration_cast<microseconds>(wtime_end-wtime_begin);
    
  printf("elapsed time = %.3lf msec\n", 0.001*elapsed_usec.count());

};

