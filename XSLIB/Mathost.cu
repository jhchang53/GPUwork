/*
	class member of Mater to transfer data to device
*/
#include <stdio.h>
#include <assert.h>

#include "Mater.h"
#include <cuda_runtime.h>
#include "Matdev.h"

void Mater::sendMat()
{
  cudaError_t rc;
  /*  nmat -> nmater  */
  rc = cudaMemcpyToSymbol(nmater,&nmat,sizeof(int));
  assert(rc == cudaSuccess);

  /*  index_iso  */
  printf("nmat=%d\n",nmat);
  for(int m=0; m <= nmat; m++) printf("  %d",ind_iso[m]);
  printf("\n");
  int *h_index_iso;
  size_t bsize_ind_iso = (nmat+1)*sizeof(int);
  rc = cudaMalloc(&h_index_iso,bsize_ind_iso);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_index_iso,ind_iso,bsize_ind_iso,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(index_iso,&h_index_iso,sizeof(int *));
  assert(rc == cudaSuccess);

  /*  iso_vec  */
  int len_iso = ind_iso[nmat];
  int *h_iso;
  size_t bsize_iso = len_iso*sizeof(int);
  rc = cudaMalloc(&h_iso,bsize_iso);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_iso,&iso_vec[0],bsize_iso,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(iso,&h_iso,sizeof(int *));
  assert(rc == cudaSuccess);

  /*  index_fmac  */
  printf("nmat=%d\n",nmat);
  for(int m=0; m <= nmat; m++) printf("  %d",ind_fmac[m]);
  printf("\n");
  int *h_index_fmac;
  size_t bsize_ind_fmac = (nmat+1)*sizeof(int);
  rc = cudaMalloc(&h_index_fmac,bsize_ind_fmac);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_index_fmac,ind_fmac,bsize_ind_fmac,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(index_fmac,&h_index_fmac,sizeof(int *));
  assert(rc == cudaSuccess);
  /*  fast_mac_vec  to fast_sigtot */
  int len_fmac = ind_fmac[nmat];
  int *h_fast_sigtot;
  size_t bsize_fmac = len_fmac*sizeof(double);
  rc = cudaMalloc(&h_fast_sigtot,bsize_fmac);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_fast_sigtot,&fast_mac_vec[0],bsize_fmac,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(fast_sigtot,&h_fast_sigtot,sizeof(double *));
  assert(rc == cudaSuccess);
  /*  index_fcdf  */
  int *h_index_fcdf;
  size_t bsize_ind_fcdf = (nmat+1)*sizeof(int);
  rc = cudaMalloc(&h_index_fcdf,bsize_ind_fcdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_index_fcdf,ind_fcdf,bsize_ind_fcdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(index_fcdf,&h_index_fcdf,sizeof(int *));
  assert(rc == cudaSuccess);
  /*  fast_cdf_vec  to fast_cdf  */
  int len_fcdf = ind_fcdf[nmat];
  int *h_fast_cdf;
  size_t bsize_fcdf = len_fcdf*sizeof(double);
  rc = cudaMalloc(&h_fast_cdf,bsize_fcdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_fast_cdf,&fast_cdf_vec[0],bsize_fcdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(fast_cdf,&h_fast_cdf,sizeof(double *));
  assert(rc == cudaSuccess);

  /*  index_tmac  */
  printf("nmat=%d\n",nmat);
  for(int m=0; m <= nmat; m++) printf("  %d",ind_tmac[m]);
  printf("\n");
  int *h_index_tmac;
  size_t bsize_ind_tmac = (nmat+1)*sizeof(int);
  rc = cudaMalloc(&h_index_tmac,bsize_ind_tmac);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_index_tmac,ind_tmac,bsize_ind_tmac,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(index_tmac,&h_index_tmac,sizeof(int *));
  assert(rc == cudaSuccess);
  /*  therm_mac_vec  to therm_sigtot */
  int len_tmac = ind_tmac[nmat];
  int *h_therm_sigtot;
  size_t bsize_tmac = len_tmac*sizeof(double);
  rc = cudaMalloc(&h_therm_sigtot,bsize_tmac);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_therm_sigtot,&therm_mac_vec[0],bsize_tmac,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(therm_sigtot,&h_therm_sigtot,sizeof(double *));
  assert(rc == cudaSuccess);
  /*  index_tcdf  */
  int *h_index_tcdf;
  size_t bsize_ind_tcdf = (nmat+1)*sizeof(int);
  rc = cudaMalloc(&h_index_tcdf,bsize_ind_tcdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_index_tcdf,ind_tcdf,bsize_ind_tcdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(index_tcdf,&h_index_tcdf,sizeof(int *));
  assert(rc == cudaSuccess);
  /*  therm_cdf_vec  to therm_cdf  */
  int len_tcdf = ind_tcdf[nmat];
  int *h_therm_cdf;
  size_t bsize_tcdf = len_tcdf*sizeof(double);
  rc = cudaMalloc(&h_therm_cdf,bsize_tcdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_therm_cdf,&therm_cdf_vec[0],bsize_tcdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(therm_cdf,&h_therm_cdf,sizeof(double *));
  assert(rc == cudaSuccess);


};	// Mater::sendMat
