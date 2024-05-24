/*
	XSLhost.cpp
*/
#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>

#include "XS.h"
#include "XSLhost.h"
#include "XSLdev.h"

void XSLhost(XSL *xsl)
{
  XSLhost_head(xsl);
  XSLhost_fast_sig(xsl);
  XSLhost_fast_xmu(xsl);
  XSLhost_fast_inel(xsl);
  XSLhost_fast_n2n(xsl);
  XSLhost_fast_n3n(xsl);
  XSLhost_therm(xsl);
  XSLhost_gamma(xsl);
  XSLhost_gprod(xsl);
  XSLhost_kerma(xsl);
};

void XSLhost_head(XSL *xsl)
{
  printf("num_fast_gps=%d num_therm_gps=%d num_gam_gps=%d num_iso=%d\n",
	xsl->num_fast_gps, xsl->num_therm_gps, xsl->num_gam_gps, xsl->num_iso);
  int fast_grps = xsl->num_fast_gps;
  int therm_grps = xsl->num_therm_gps;
  int neut_grps = xsl->num_fast_gps + xsl->num_therm_gps;
  int gam_grps = xsl->num_gam_gps;
  int niso = xsl->num_iso;

  cudaError_t rc;
  
  rc = cudaMemcpyToSymbol(num_fast_gps,&fast_grps,sizeof(int));
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(num_therm_gps,&therm_grps,sizeof(int));
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(num_neut_gps,&neut_grps,sizeof(int));
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(num_gam_gps,&gam_grps,sizeof(int));
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(num_iso,&niso,sizeof(int));
  assert(rc == cudaSuccess);

  double *h_ev_neut;
  size_t bsize_neut = (neut_grps+1)*sizeof(double);
  rc = cudaMalloc(&h_ev_neut,bsize_neut);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ev_neut,xsl->ev_neut,bsize_neut,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ev_neut,&h_ev_neut,sizeof(double *));
  assert(rc == cudaSuccess);

  double *h_ev_gam;
  size_t bsize_gam = (gam_grps+1)*sizeof(double);
  rc = cudaMalloc(&h_ev_gam,bsize_gam);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ev_gam,xsl->ev_gam,bsize_gam,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ev_gam,&h_ev_gam,sizeof(double *));
  assert(rc == cudaSuccess);

  int *h_id_fast;
  size_t bsize_id_fast = niso*sizeof(int);
  rc = cudaMalloc(&h_id_fast,bsize_id_fast);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_id_fast,xsl->id_fast,bsize_id_fast,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(id_fast,&h_id_fast,sizeof(int *));
  assert(rc == cudaSuccess);

  int *h_id_therm;
  size_t bsize_id_therm = niso*sizeof(int);
  rc = cudaMalloc(&h_id_therm,bsize_id_therm);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_id_therm,xsl->id_therm,bsize_id_therm,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(id_therm,&h_id_therm,sizeof(int *));
  assert(rc == cudaSuccess);

  int *h_id_gam;
  size_t bsize_id_gam = niso*sizeof(int);
  rc = cudaMalloc(&h_id_gam,bsize_id_gam);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_id_gam,xsl->id_gam,bsize_id_gam,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(id_gam,&h_id_gam,sizeof(int *));
  assert(rc == cudaSuccess);
};	// XSLhost_head

int XSLhost_fast_sig(XSL *xsl)
{
  /*  makee cdf of fast reaction type  */
  int niso = xsl->num_iso;
  double *amass = new double[niso];
  int ngrps = xsl->num_fast_gps;
  int indmax = niso*ngrps*5;	// prepare abs, el, inel, n2n, n3n
  int *indfast = new int[niso+1];
  int look = -1; // niso-1;
  int ind = 0;
  double *fastcdf = new double[indmax];
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    amass[iso] = fast.awr;	// collect atomic mass
    /*  create  fast reaction type cdf   */
    double cdf[500];	// reserve max fast
    for(int g=0; g < ngrps; g++) {
      cdf[5*g] = fast.xsabs[g]/fast.sigtot[g];
    }
    for(int g=0; g < ngrps; g++) {
      cdf[5*g+1] = cdf[5*g]+fast.elas[g]/fast.sigtot[g];
    }
    int inel_grps = fast.inel_grps;
    for(int g=0; g < inel_grps; g++) {
      cdf[5*g+2] = cdf[5*g+1]+fast.sig_inel[g]/fast.sigtot[g];
    }
    for(int g=inel_grps; g < ngrps; g++) {
      cdf[5*g+2] = cdf[5*g+1];
    }
    int n2n_grps = fast.n2n_grps;
    for(int g=0; g < n2n_grps; g++) {
      cdf[5*g+3] = cdf[5*g+2]+fast.sig_n2n[g]/fast.sigtot[g];
    }
    for(int g=n2n_grps; g < ngrps; g++) {
      cdf[5*g+3] = cdf[5*g+2];
    }
    /*  enforce last cdf for n3n as 1.0  */
    for(int g=0; g < ngrps; g++) {
      cdf[5*g+4] = 1.0;
    }
    if(iso==look) {
      printf("cdf");
      for(int g=0; g < ngrps; g++) {
        for(int i=0; i < 5; i++) printf(" %.12lf",cdf[5*g+i]);
        printf("\n");
      }
    }
    /*  create full  size is 5*ngrp  */
    indfast[iso] = ind;
    for(int g5=0; g5 < 5*ngrps; g5++) {
      fastcdf[ind+g5] = cdf[g5];
    }
    ind += 5*ngrps;
  }	// for iso
  indfast[niso] = ind;
  /*  send to GPU  */
  cudaError_t rc;

  double *h_awr;
  size_t bsize_awr = niso*sizeof(double);
  rc = cudaMalloc(&h_awr,bsize_awr);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_awr,amass,bsize_awr,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(awr,&h_awr,sizeof(double *));
  assert(rc == cudaSuccess);

  int *h_ind_fast;
  size_t bsize_ind = (niso+1)*sizeof(int);
  rc = cudaMalloc(&h_ind_fast,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_fast,indfast,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_fast,&h_ind_fast,sizeof(int *));
  assert(rc == cudaSuccess);

  double *h_fast_cdf;
  size_t bsize_cdf = niso*5*ngrps*sizeof(double);
  rc = cudaMalloc(&h_fast_cdf,bsize_cdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_fast_cdf,fastcdf,bsize_cdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(fast_cdf,&h_fast_cdf,sizeof(double *));
  assert(rc == cudaSuccess);
  return ind;
};	//  XSLhost_fast_sig

