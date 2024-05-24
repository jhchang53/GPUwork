/*
	XSLhost2.cu - part 2
*/
#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>

#include "XS.h"
#include "XSLhost.h"
#include "XSLdev.h"

int XSLhost_fast_xmu(XSL *xsl)
{
  /*  elastic xmu  */
  int niso = xsl->num_iso;
  int nfast_gps = xsl->num_fast_gps;
  int *indxmu = new int[niso*nfast_gps+1];
  int ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    for(int g=0; g < nfast_gps; g++) {
      indxmu[iso*nfast_gps+g] = ind;
      ind += fast.nxmus[g];
    }
  }
  int indlast = ind;
  indxmu[niso*nfast_gps] = indlast;
  int look=-1; // niso-1;
  if(look >=0 ) {
    printf(" indmus:");
    for(int g=0; g < nfast_gps; g++) printf(" %d",indxmu[look*nfast_gps+g]);
    printf("\n");
  }
  /*  collect xmu  */
  double *xmus = new double[indlast];
  ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    for(int g=0; g < nfast_gps; g++) {
      assert(indxmu[iso*nfast_gps+g] == ind);
      int nxmu = fast.nxmus[g];
      for(int i=0; i < nxmu; i++) {
        xmus[ind+i] = fast.xmu[g][i];
      }
      ind += fast.nxmus[g];
    }
  }
  if(look >= 0) {
    for(int g=0; g < nfast_gps; g++) {
      printf("xmu%d:",g);
      int idx = indxmu[look*nfast_gps+g];
      int ndx = indxmu[look*nfast_gps+g+1];
      for(int i=idx; i < ndx; i++) {
        printf(" %.4lf",xmus[i]);
      }
      printf("\n");
    }
  }
  /*  send to GPU  */
    cudaError_t rc;
  int *h_ind_xmu;
  size_t bsize_ind = (niso*nfast_gps+1)*sizeof(int);
  rc = cudaMalloc(&h_ind_xmu,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_xmu,indxmu,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_xmu,&h_ind_xmu,sizeof(int *));
  assert(rc == cudaSuccess);

  int *h_xmu;
  size_t bsize_xmu = indlast*sizeof(double);
  rc = cudaMalloc(&h_xmu,bsize_xmu);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_xmu,xmus,bsize_xmu,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(xmu,&h_xmu,sizeof(double *));
  assert(rc == cudaSuccess);
  return indlast;
};	// XSLhost_fast_xmu

int XSLhost_fast_inel(XSL *xsl)
{
  /*  collect inelastic  cdf and aniso . both have size inel_cdfs  */
  int niso = xsl->num_iso;
  int *indinel = new int[niso];
  int *inelgrps = new int[niso];
  int *inelcdfs = new int[niso];
  int ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    int grps = fast.inel_grps;
    int cdfs = fast.inel_cdfs;
    inelgrps[iso] = grps;
    inelcdfs[iso] = cdfs;
    indinel[iso] = ind;
    ind += grps*(cdfs+1);
  }
  int indlast = ind;
  double *cdf = new double[indlast];
  double *aniso = new double[indlast];
  ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    int grps = fast.inel_grps;
    int cdfs = fast.inel_cdfs;
    int lcdf = cdfs+1;
    assert(ind == indinel[iso]);
    for(int ig=0; ig < grps; ig++) {
      for(int jg=0; jg < lcdf; jg++) {
        cdf[ind+ig*lcdf+jg] = fast.cdf_inel[ig*lcdf+jg];
        if(fast.aniso_inel == nullptr) aniso[ind+ig*lcdf+jg] = 0.0;
        else aniso[ind+ig*lcdf+jg] = fast.aniso_inel[ig*lcdf+jg];
      }
    }	// for ig
    ind += grps*(cdfs+1);
  }	// for iso
  /*  send to GPU  */
  cudaError_t rc;
  int *h_ind_inel;
  size_t bsize_ind = niso*sizeof(int);
  rc = cudaMalloc(&h_ind_inel,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_inel,indinel,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_inel,&h_ind_inel,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_inel_grps;
  rc = cudaMalloc(&h_inel_grps,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_inel_grps,inelgrps,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(inel_grps,&h_inel_grps,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_inel_cdfs;
  rc = cudaMalloc(&h_inel_cdfs,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_inel_cdfs,inelcdfs,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(inel_cdfs,&h_inel_cdfs,sizeof(int *));
  assert(rc == cudaSuccess);

  int *h_cdf_inel;
  size_t bsize_cdf = indlast*sizeof(double);
  rc = cudaMalloc(&h_cdf_inel,bsize_cdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_cdf_inel,cdf,bsize_cdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(cdf_inel,&h_cdf_inel,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_aniso_inel;
  rc = cudaMalloc(&h_aniso_inel,bsize_cdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_aniso_inel,aniso,bsize_cdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(aniso_inel,&h_aniso_inel,sizeof(int *));
  assert(rc == cudaSuccess);

  return indlast;
};	// XSLhost_fast_inel

int XSLhost_fast_n2n(XSL *xsl)
{
  /*  n2n cdfs  */
  int niso = xsl->num_iso;
  int *indn2n  = new int[niso];
  int *n2ngrps = new int[niso];
  int *n2ncdfs = new int[niso];
  int ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    int grps = fast.n2n_grps;
    int cdfs = fast.n2n_cdfs;
    n2ngrps[iso] = grps;
    n2ncdfs[iso] = cdfs;
    indn2n[iso] = ind;
    ind += grps*(cdfs+1);
  }
  int indlast = ind;
  double *cdf = new double[indlast];
  ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    int grps = fast.n2n_grps;
    int cdfs = fast.n2n_cdfs;
    int lcdf = cdfs+1;
    assert(ind == indn2n[iso]);
    for(int ig=0; ig < grps; ig++) {
      for(int jg=0; jg < lcdf; jg++) {
        cdf[ind+ig*lcdf+jg] = fast.cdf_n2n[ig*lcdf+jg];
      }
    }   // for ig
    ind += grps*(cdfs+1);
  }     // for iso
  /*  send to GPU  */
  cudaError_t rc;
  int *h_ind_n2n;
  size_t bsize_ind = niso*sizeof(int);
  rc = cudaMalloc(&h_ind_n2n,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_n2n,indn2n,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_n2n,&h_ind_n2n,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_n2n_grps;
  rc = cudaMalloc(&h_n2n_grps,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_n2n_grps,n2ngrps,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(n2n_grps,&h_n2n_grps,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_n2n_cdfs;
  rc = cudaMalloc(&h_n2n_cdfs,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_n2n_cdfs,n2ncdfs,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(n2n_cdfs,&h_n2n_cdfs,sizeof(int *));
  assert(rc == cudaSuccess);

  int *h_cdf_n2n;
  size_t bsize_cdf = indlast*sizeof(double);
  rc = cudaMalloc(&h_cdf_n2n,bsize_cdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_cdf_n2n,cdf,bsize_cdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(cdf_n2n,&h_cdf_n2n,sizeof(int *));
  assert(rc == cudaSuccess);
  return indlast;
};	// XSLhost_fast_n2n

int XSLhost_fast_n3n(XSL *xsl)
{
  /*  n3n cdfs  */
  int niso = xsl->num_iso;
  int *indn3n  = new int[niso];
  int *n3ngrps = new int[niso];
  int *n3ncdfs = new int[niso];
  int ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    int grps = fast.n3n_grps;
    int cdfs = fast.n3n_cdfs;
    n3ngrps[iso] = grps;
    n3ncdfs[iso] = cdfs;
    indn3n[iso] = ind;
    ind += grps*(cdfs+1);
  }
  int indlast = ind;
  double *cdf = new double[indlast];
  ind = 0;
  for(int iso=0; iso < niso; iso++) {
    XSfast fast = xsl->xs[iso].fast;
    int grps = fast.n3n_grps;
    int cdfs = fast.n3n_cdfs;
    int lcdf = cdfs+1;
    assert(ind == indn3n[iso]);
    for(int ig=0; ig < grps; ig++) {
      for(int jg=0; jg < lcdf; jg++) {
        cdf[ind+ig*lcdf+jg] = fast.cdf_n3n[ig*lcdf+jg];
      }
    }   // for ig
    ind += grps*(cdfs+1);
  }     // for iso
  /*  send to GPU  */
  cudaError_t rc;
  int *h_ind_n3n;
  size_t bsize_ind = niso*sizeof(int);
  rc = cudaMalloc(&h_ind_n3n,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_ind_n3n,indn3n,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(ind_n3n,&h_ind_n3n,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_n3n_grps;
  rc = cudaMalloc(&h_n3n_grps,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_n3n_grps,n3ngrps,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(n3n_grps,&h_n3n_grps,sizeof(int *));
  assert(rc == cudaSuccess);
  int *h_n3n_cdfs;
  rc = cudaMalloc(&h_n3n_cdfs,bsize_ind);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_n3n_cdfs,n3ncdfs,bsize_ind,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(n3n_cdfs,&h_n3n_cdfs,sizeof(int *));
  assert(rc == cudaSuccess);

  int *h_cdf_n3n;
  size_t bsize_cdf = indlast*sizeof(double);
  rc = cudaMalloc(&h_cdf_n3n,bsize_cdf);
  assert(rc == cudaSuccess);
  rc = cudaMemcpy(h_cdf_n3n,cdf,bsize_cdf,cudaMemcpyHostToDevice);
  assert(rc == cudaSuccess);
  rc = cudaMemcpyToSymbol(cdf_n3n,&h_cdf_n3n,sizeof(int *));
  assert(rc == cudaSuccess);
  return indlast;
};	// XSLhost_fast_n3n

