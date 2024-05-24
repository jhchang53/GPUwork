/*
	xslib.cu
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>

#include "XSLIB.h"
#include "XSLhost.h"
#include "XSLdev.h"

int main()
{
  XSLIB *xslib = new XSLIB();
  XSL *xsl = xslib->openlib("seraMC.libupd4x.punch");
  XSLhost(xsl);

  /*  call GPU  */
  XSLset<<<1,1>>>(8);
  
  cudaDeviceSynchronize();
  delete xslib;
};
