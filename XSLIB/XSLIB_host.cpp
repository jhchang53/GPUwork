/*
	XSLIBhost.cpp
*/
#include <stdio.h>
#include <assert.h>
#include "XSLIB.h"
#include <cuda_runtime.h>

__device__ XS *d_xs;

void XSLIB::sendToDevice(XS *xs)
{
  cudaError_t rc;
  const XS *h_xs;
  rc = cudaMalloc((void **)&h_xs,sizeof(XS));
  assert(rc == cudaSuccess);
  XS*  ptr;
  rc = cudaGetSymbolAddress((void **)&ptr, h_xs);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(ptr,xs,sizeof(xs),cudaMemcpyHostToDevice);

  // rc = cudaMemcpy(h_xs,xs,sizeof(XS), cudaMemcpyHostToDevice);
  // assert(rc == cudaSuccess);
  // rc = cudaMemcpyToSymbol(d_xs,&h_xs,sizeof(XS *),0,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);

};	// XSLIB::sendToDevice
