/*
	XSLhost3.cu
*/
#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>

#include "XS.h"
#include "XSLhost.h"
#include "XSLdev.h"

int XSLhost_therm(XSL *xsl)
{
  /*  make thermal reaction type cdf  */
  /*  only capture or scattering  */
  int niso = xsl->num_iso;
  int ngrps = xsl->num_therm_gps;
  int *th_jscat = new int[niso];
  int *th_nbound = new int[niso];
  /*  we dont need index but... */
  int *indtherm = new int[niso+1];
  int ind  = 0;
  double *thermcdf = new double[niso*ngrps*2];
  for(int iso=0; iso < niso; iso++) {
    XStherm therm = xsl->xs[iso].therm;
    th_jscat[iso]  = therm.jscat;
    th_nbound[iso] = therm.nbound;
    /*  thermal can be ommitted for Hydrogen  508 0  */
    if(therm.sigtot == nullptr) {
      for(int g=0; g < ngrps; g++) {
        thermcdf[ind+g*2] = 0.0;	// all is scattering
        thermcdf[ind+g*2+1] = 1.0;
      }
    }
    else {
      for(int g=0; g < ngrps; g++) {
        double capture = therm.siga[g]/therm.sigtot[g];
        thermcdf[ind+g*2] = capture;
        thermcdf[ind+g*2+1] = 1.0;	// enforce
      }
    }
    indtherm[iso] = ind;
    ind += 2*ngrps;
  }
  int indlast = ind;
  indtherm[niso] = indlast;
  /*  send to GPU  */
  cudaError_t rc;

  int *h_jscat;
  size_t bsize_jscat = niso*sizeof(int);
  rc = cudaMalloc(&h_jscat,bsize_jscat);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_jscat,th_jscat,bsize_jscat,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(jscat,&h_jscat,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_nbound;
  size_t bsize_nbound = niso*sizeof(int);
  rc = cudaMalloc(&h_nbound,bsize_nbound);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_nbound,th_nbound,bsize_nbound,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(nbound,&h_nbound,sizeof(int *));
  assert(rc == cudaSuccess);

  int *h_ind_therm;
  size_t bsize_ind = (niso+1)*sizeof(int);
  rc = cudaMalloc(&h_ind_therm,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_therm,indtherm,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_therm,&h_ind_therm,sizeof(int *));
  assert(rc == cudaSuccess);

  double *h_therm_cdf;
  size_t bsize_cdf = niso*2*ngrps*sizeof(double);
  rc = cudaMalloc(&h_therm_cdf,bsize_cdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_therm_cdf,thermcdf,bsize_cdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(therm_cdf,&h_therm_cdf,sizeof(double *));
  assert(rc == cudaSuccess);

  return indlast;
};	// XSLhost_therm

int XSLhost_gprod(XSL *xsl)
{
  /*  capture gamma  */
  int niso = xsl->num_iso;
  double *sig2200  = new double[niso];
  double *siggprod = new double[niso];
  int *numline = new int[niso];
  int *indgprod = new int[niso+1];
  int ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSgprod gprod = xsl->xs[iso].gprod;
    sig2200[iso] = gprod.sig_2200;
    siggprod[iso] = gprod.sig_gprod;
    numline[iso] = gprod.num_gamline;
    indgprod[iso] = ind;
    ind += numline[iso];
  }
  int indlast = ind;
  indgprod[niso] = indlast;
  double *yield = new double[indlast];
  double *ener  = new double[indlast];
  ind  = 0;
  for(int iso=0; iso < niso; iso++) {
    XSgprod gprod = xsl->xs[iso].gprod;
    int numl = numline[iso];
    for(int n=0; n < numl; n++) {
      yield[ind+n] = gprod.g_yield[n];
      ener[ind+n]  = gprod.g_energy[n];
    }
    ind += numl;
  }
  /*  send to GPU  */
  cudaError_t rc;

  double *h_sig_2200;
  size_t bsize_sig = niso*sizeof(double);
  rc = cudaMalloc(&h_sig_2200,bsize_sig);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_sig_2200,sig2200,bsize_sig,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(sig_2200,&h_sig_2200,sizeof(double *));
  assert(rc == cudaSuccess);
  double *h_sig_gprod;
  rc = cudaMalloc(&h_sig_gprod,bsize_sig);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_sig_gprod,siggprod,bsize_sig,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(sig_gprod,&h_sig_gprod,sizeof(double *));
  assert(rc == cudaSuccess);
  int *h_num_gamline;
  size_t bsize_numl = niso*sizeof(int);
  rc = cudaMalloc(&h_num_gamline,bsize_numl);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_num_gamline,numline,bsize_numl,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(num_gamline,&h_num_gamline,sizeof(double *));
  assert(rc == cudaSuccess);

  int *h_ind_gprod;
  size_t bsize_ind = (niso+1)*sizeof(int);
  rc = cudaMalloc(&h_ind_gprod,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_gprod,indgprod,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_gprod,&h_ind_gprod,sizeof(int *));
  assert(rc == cudaSuccess);

  double *h_g_yield;
  size_t bsize_yield = indlast*sizeof(double);
  rc = cudaMalloc(&h_g_yield,bsize_yield);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_g_yield,yield,bsize_yield,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(g_yield,&h_g_yield,sizeof(double *));
  assert(rc == cudaSuccess);

  double *h_g_energy;
  rc = cudaMalloc(&h_g_energy,bsize_yield);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_g_energy,ener,bsize_yield,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(g_energy,&h_g_energy,sizeof(double *));
  assert(rc == cudaSuccess);

  return indlast;
};      // XSLhost_gprod

