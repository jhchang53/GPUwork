/*
	test driver

*/
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <curand.h>

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
printf("Error at %s:%d\n",__FILE__,__LINE__); \
return EXIT_FAILURE;}} while(0)

#define THREADS_PER_BLOCK 3
#define BLOCK_COUNT 4

__global__ void setup_kernel(curandStateMRG32k3a *state)
{
  int id = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  unsigned long long seed = id+293;
  unsigned long long sequence = 1;
  unsigned long long offset = 71;
  curand_init(seed,sequence,offset,&state[id]);
};

__global__ void generate_kernel(curandStateMRG32k3a *state)
{
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  // unsigned int x;
  double x;
  /* Copy state to local memory for efficiency */
  curandStateMRG32k3a localState = state[id];
  /* Generate pseudo-random unsigned ints */
  for(int i = 0; i < 4; i++) {
    // x = curand(&localState);
    // printf(" id=%d x=%d\n",id,x);
    x = curand_uniform_double(&localState);
    printf(" id=%d x=%.5lf\n",id,x);
  }
  /* Copy state back to global memory */
  state[id] = localState;
};


int main()
{
  const unsigned int threadsPerBlock = THREADS_PER_BLOCK;
  const unsigned int blockCount = BLOCK_COUNT;
  const unsigned int totalThreads = threadsPerBlock * blockCount;
  /*  initialize */
  curandStateMRG32k3a *devMRGStates;
  CUDA_CALL(cudaMalloc((void **)&devMRGStates, totalThreads *
	sizeof(curandStateMRG32k3a)));
  setup_kernel<<<BLOCK_COUNT,THREADS_PER_BLOCK>>>(devMRGStates);
  for(int iter=0; iter < 2; iter++) {
  generate_kernel<<<BLOCK_COUNT,THREADS_PER_BLOCK>>>(devMRGStates);
  cudaDeviceSynchronize();
    printf(" iter=%d done\n",iter);
  }
};

