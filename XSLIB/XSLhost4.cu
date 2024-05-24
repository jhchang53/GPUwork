/*
	XSLhost4. cu
*/
#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>

#include "XS.h"
#include "XSLhost.h"
#include "XSLdev.h"

int XSLhost_gamma(XSL *xsl)
{
  /*  make gamma reaction type cdf  */
  /*  photoabsorption, Compton, pair production  */
  int niso = xsl->num_iso;
  int ngrps = xsl->num_gam_gps;
  /*  we dont need index but... */
  int *indgam = new int[niso+1];
  int ind  = 0;
  double *gamcdf = new double[niso*ngrps*3];
  for(int iso=0; iso < niso; iso++) {
    XSgamma gam = xsl->xs[iso].gamma;
    for(int g=0; g < ngrps; g++) {
      gamcdf[ind+g*3] = gam.phab[g]/gam.sigtot[g];
      gamcdf[ind+g*3+1] = (gam.phab[g]+gam.compt[g])/gam.sigtot[g];
      gamcdf[ind+g*3+2] = 1.0;      // enforce
    }
    indgam[iso] = ind;
    ind += 3*ngrps;
  }
  int indlast = ind;
  indgam[niso] = indlast;
  /*  send to GPU  */
  cudaError_t rc;
  int *h_ind_gam;
  size_t bsize_ind = (niso+1)*sizeof(int);
  rc = cudaMalloc(&h_ind_gam,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_gam,indgam,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_gam,&h_ind_gam,sizeof(int *));
  assert(rc == cudaSuccess);

  double *h_gam_cdf;
  size_t bsize_cdf = niso*3*ngrps*sizeof(double);
  rc = cudaMalloc(&h_gam_cdf,bsize_cdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_gam_cdf,gamcdf,bsize_cdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(gam_cdf,&h_gam_cdf,sizeof(double *));
  assert(rc == cudaSuccess);

  return indlast;
};	// XSLhost_gamma

int XSLhost_kerma(XSL *xsl)
{
  /*  kerma  */
  int niso = xsl->num_iso;
  int ngrps_n = xsl->num_fast_gps + xsl->num_therm_gps;
  double *kerma_n = new double[niso*ngrps_n];
  int ngrps_g = xsl->num_gam_gps;
  double *kerma_g = new double[niso*ngrps_g];
  
  for(int iso=0; iso < niso; iso++) {
    XSkerma kerma = xsl->xs[iso].kerma;
    int ind_n = iso*ngrps_n;
    for(int g=0; g < ngrps_n; g++) kerma_n[ind_n+g] = kerma.xkerma_n[g];
    int ind_g = iso*ngrps_g;
    for(int g=0; g < ngrps_g; g++) kerma_g[ind_g+g] = kerma.xkerma_g[g];
  }
  /*  send to GPU  */
  cudaError_t rc;

  double *h_xkerma_n;
  size_t bsize_n = niso*ngrps_n*sizeof(double);
  rc = cudaMalloc(&h_xkerma_n,bsize_n);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_xkerma_n,kerma_n,bsize_n,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(xkerma_n,&h_xkerma_n,sizeof(double *));

  double *h_xkerma_g;
  size_t bsize_g = niso*ngrps_g*sizeof(double);
  rc = cudaMalloc(&h_xkerma_g,bsize_g);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_xkerma_g,kerma_g,bsize_g,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(xkerma_g,&h_xkerma_g,sizeof(double *));

  return niso*(ngrps_n+ngrps_g);
};	// XSLhost_kerma
