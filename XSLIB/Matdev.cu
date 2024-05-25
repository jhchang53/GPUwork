/*
	Matdev.cu
*/
#include <stdio.h>
#include <assert.h>
#define __Matdev
#include "Matdev.h"

__global__ void MatSet()
{
  int prlev = 0;
  printf(" nmater=%d\n",nmater);
  printf("  index \n");
  for(int m=0; m <= nmater; m++) {
    printf(" %d",index_iso[m]);
  }
  printf("\n");
  printf(" iso:");
  for(int m=0; m < nmater; m++) {
    for(int i=index_iso[m]; i < index_iso[m+1]; i++) {
      printf(" %d",iso[i]);
    }
    printf("\n");
  }
  int nfgrps = 10;
  for(int m=0; m < nmater; m++) {
    printf("fast%d:",m);
    for(int i=index_fmac[m]; i < index_fmac[m+1]; i++) {
      printf(" %.5le",fast_sigtot[i]);
      if(i%6==5) printf("\n  ");
    }
    printf("\n");
  }
  for(int m=0; m < nmater; m++) {
    printf("fcdf%d:",m);
    int niso =  index_iso[m+1]-index_iso[m];
    for(int i=index_iso[m]; i < index_iso[m+1]; i++) {
      printf(" %5d",iso[i]);
    }
    printf("\n");
    int ind = index_fcdf[m];
    if(index_fcdf[m+1] > ind) {
      for(int g=0; g < nfgrps; g++) {
        printf("cdf%d:",g);
        for(int j=0; j < niso; j++) {
          printf(" %.5lf",fast_cdf[ind+g*niso+j]);
        }
        printf("\n");
      }
    }
  }
  int ntgrps = 10;
  for(int m=0; m < nmater; m++) {
    printf("therm%d:",m);
    for(int i=index_tmac[m]; i < index_tmac[m+1]; i++) {
      printf(" %.5le",therm_sigtot[i]);
      if(i%6==5) printf("\n  ");
    }
    printf("\n");
  }
  for(int m=0; m < nmater; m++) {
    printf("tcdf%d:",m);
    int niso =  index_iso[m+1]-index_iso[m];
    for(int i=index_iso[m]; i < index_iso[m+1]; i++) {
      printf(" %5d",iso[i]);
    }
    printf("\n");
    int ind = index_tcdf[m];
    if(index_tcdf[m+1] > ind) {
      for(int g=0; g < ntgrps; g++) {
        printf("cdf%d:",g);
        for(int j=0; j < niso; j++) {
          printf(" %.5lf",therm_cdf[ind+g*niso+j]);
        }
        printf("\n");
      }
    }
  }

};
