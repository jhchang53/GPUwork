/*
	B10kerma.cu
*/
#include <stdio.h>

__device__ int b10_size;
__device__ double *b10_ener,*b10_kerm;

__global__ void B10kerma(double En)
{
  /*  Wheeler algorithm for increasing ordered data  */
  /*  we take linear extrapolation at both ends  */
  int ifrom = 0;
  int ito = b10_size-1;
  while(ito-ifrom > 1) {
    int iguess = ifrom + (ito-ifrom)/2;
    if(En < b10_ener[iguess]) ito = iguess;
    else ifrom = iguess;
  }
  int ig = ifrom;
  printf("E=%.2le size=%d ig=%d E0=%.2le %.2le\n",
    En,b10_size,ig,b10_ener[ig],b10_ener[ig+1]);
  double weit = (b10_ener[ig+1]-En)/(b10_ener[ig+1]-b10_ener[ig]);
  double kerma = weit*b10_kerm[ig] + (1.0-weit)*b10_kerm[ig+1];
  printf(" weit=%.4lf  kerma=%.5le\n",weit,kerma);
};
