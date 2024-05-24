/*
	track.h
*/
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <curand.h>
#define THREADS_PER_BLOCK 3
#define BLOCK_COUNT 2

__global__ void setup_kernel(curandStateMRG32k3a *state,
  int ncgs, int *d_CGtype, double *d_CGparm,
  int *d_nxyz, double *d_dxyz, char *d_uv,  float *d_tally,
  double Rsrc, int nsrcg, double *d_SRCen, double *d_SRCcdf,
  double target[], double Zb, double theta, double phi, double psi);


__global__ void generate_kernel();



