#include <stdio.h>
#include "Xor128.h"

int main()
{
  init_xor128();
  for(int j=0; j < 3; j++) jump_xor128();
  int nran = 100000;
  double sum1 = 0;      double sum11 = 0;
  double sum2 = 0;      double sum22 = 0;
  double sum12 = 0;
  for(int n=0; n < nran; n++) {
    double ran1 = ranf();
    double ran2 = ranf();
    sum1 += ran1;
    sum11 += (ran1-0.5)*(ran1-0.5);
    sum2 += ran2;
    sum22 += (ran2-0.5)*(ran2-0.5);
    sum12 += (ran1-0.5)*(ran2-0.5);
  }
  double var11 = sum11/(nran-1);
  double var22 = sum22/(nran-1);
  double var12 = sum12/(nran-1);
  printf(" avg1=%.5lf var11=%.4lf\n",sum1/nran,12*var11);
  printf(" avg1=%.5lf var22=%.4lf\n",sum2/nran,12*var22);
  printf(" var12=%.4lf\n",12*var12);
};

