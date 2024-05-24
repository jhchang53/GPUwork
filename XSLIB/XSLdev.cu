/*
  XSLset.cu
*/
#include <stdio.h>
#include <assert.h>
#define __XSLdev
#include "XSLdev.h"

__global__ void XSLset(int look)
{
  int prlev = 0;
  printf("=== XSLset at Device ===\n");
  printf(" fast=%d therm=%d neut=%d gam=%d iso=%d\n",
	num_fast_gps,num_therm_gps,num_neut_gps,num_gam_gps,num_iso);
  if(prlev > 1) {
    printf("XSLset ev_neut:");
    for(int g=0; g <= num_neut_gps; g++) {
      printf(" %.5le",ev_neut[g]);
      if(g%6==5) printf("\n  ");
    }
    printf("\n");

    printf("XSLset ev_gam:");
    for(int g=0; g <= num_gam_gps; g++) {
      printf(" %.5le",ev_gam[g]);
      if(g%6==5) printf("\n  ");
    }
    printf("\n");
  }
  if(prlev) {
    printf("== id_fast therm gam  jscat  nbound\n");
    for(int iso=0; iso < num_iso; iso++) {
      printf("   %5d %5d %5d  %5d %5d\n",id_fast[iso],id_therm[iso],id_gam[iso],
	jscat[iso],nbound[iso]);
    }
  }
  printf(" look=%d awr=%.5lf\n",look,awr[look]);

  // checkFast(look);
  // checkXmu(look);
  // checkInel(look);
  // checkN2n(look);
  // checkN3n(look);
  // checkTherm(look);
  // checkGam(look);
  // checkGprod(look);
  // checkKerma(look);
};	// XSLset

__device__ void checkFast(int iso)
{
  int indx = ind_fast[iso];
  int ngrps = num_fast_gps;
  printf("fast_cdf:");
  for(int g=0; g < ngrps; g++) {
    for(int i=0; i < 5; i++) printf(" %.6lf",fast_cdf[indx+5*g+i]);
    printf("\n");
  }
};	//  checkFast

__device__ void checkXmu(int iso)
{
  printf("  checkXmu %d\n",iso);
  int ngrps = num_fast_gps;
  for(int g=0; g < ngrps; g++) {
    printf("xmu%d:",g);
    int idx = ind_xmu[iso*ngrps+g];
    int ndx = ind_xmu[iso*ngrps+g+1];
    for(int i=idx; i < ndx; i++) {
      printf(" %.4lf",xmu[i]);
    }
    printf("\n");
  }
};	//  checkXmu

__device__ void checkInel(int iso)
{
  printf("  checkInel %d\n",iso);
  int lcdf = inel_cdfs[iso]+1;
  int ind = ind_inel[iso];
  for(int ig=0; ig < inel_grps[iso]; ig++) {
    printf("cdf%d:",ig);
    for(int jg=0; jg < lcdf; jg++) {
      printf(" %.5le",cdf_inel[ind+ig*lcdf+jg]);
      if(jg%6 == 5) printf("\n  ");
    }
    printf("\n");
  }
  for(int ig=0; ig < inel_grps[iso]; ig++) {
    printf("ani%d:",ig);
    for(int jg=0; jg < lcdf; jg++) {
      printf(" %.5le",aniso_inel[ind+ig*lcdf+jg]);
      if(jg%6 == 5) printf("\n  ");
    }
    printf("\n");
  }
};	// checkInel

__device__ void checkN2n(int iso)
{
  printf("  checkN2n %d\n",iso);
  int lcdf = n2n_cdfs[iso]+1;
  int ind = ind_n2n[iso];
  for(int ig=0; ig < n2n_grps[iso]; ig++) {
    printf("cdf%d:",ig);
    for(int jg=0; jg < lcdf; jg++) {
      printf(" %.5le",cdf_n2n[ind+ig*lcdf+jg]);
      if(jg%6 == 5) printf("\n  ");
    }
    printf("\n");
  }
};      // checkN2n

__device__ void checkN3n(int iso)
{
  printf("  checkN3n %d\n",iso);
  int lcdf = n3n_cdfs[iso]+1;
  int ind = ind_n3n[iso];
  for(int ig=0; ig < n3n_grps[iso]; ig++) {
    printf("cdf%d:",ig);
    for(int jg=0; jg < lcdf; jg++) {
      printf(" %.5le",cdf_n3n[ind+ig*lcdf+jg]);
      if(jg%6 == 5) printf("\n  ");
    }
    printf("\n");
  }
};      // checkN3n

__device__ void checkTherm(int iso)
{
  printf("  checkTherm %d\n",iso);
  int lcdf = 2;
  int ind = ind_therm[iso];
  for(int g=0; g < num_therm_gps; g++) {
    printf("cdf%d:",g);
    for(int i=0; i < 2; i++) {
      printf(" %.5le",therm_cdf[ind+g*lcdf+i]);
    }
    printf("\n");
  }
};      // checkTherm

__device__ void checkGam(int iso)
{
  printf("  checkGamm %d  gam_gps=%d\n",iso,num_gam_gps);
  int lcdf = 3;
  int ind = ind_gam[iso];
  for(int g=0; g < num_gam_gps; g++) {
    printf("cdf%d:",g);
    for(int i=0; i < 3; i++) {
      printf(" %.5le",gam_cdf[ind+g*lcdf+i]);
    }
    printf("\n");
  }
};      // checkTherm

__device__ void checkGprod(int iso)
{
  printf("  checkGprod %d\n",iso);
  printf(" sig2200=%.5le  siggprod=%.5le\n",sig_2200[iso],sig_gprod[iso]);
  int ind = ind_gprod[iso];
  int numl = num_gamline[iso];
  printf("ind=%d numl=%d\n",ind,numl);
  assert(numl < 999);
  printf("   g_energy     g_yield\n");
  for(int n=0; n < numl; n++) {
    printf(" %.5le  %.5le\n",g_energy[ind+n],g_yield[ind+n]);
  }
};	// checkGprod

__device__ void checkKerma(int iso)
{
  printf("  checkKerma %d\n",iso);
  /*  kerma array is fixed size  */
  int ind_n = iso*num_neut_gps;
  printf("kerma_n:");
  for(int g=0; g < num_neut_gps; g++) {
    printf(" %.5le",xkerma_n[ind_n+g]);
    if(g%6 == 5) printf("\n   ");
  }
  printf("\n");
  int ind_g = iso*num_gam_gps;
  printf("kerma_g:");
  for(int g=0; g < num_gam_gps; g++) {
    printf(" %.5le",xkerma_g[ind_g+g]);
    if(g%6 == 5) printf("\n   ");
  }
  printf("\n");
};	// checkKerma

