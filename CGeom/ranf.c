#include <stdio.h>
#include "Xor128.h"

int main()
{
  init_xor128();
  xor128_seed_state = 0;
  xor128_s[0] = xor128_next_int();
  xor128_s[1] = xor128_next_int();
  uint64_t s1[2],s2[2];
  s1[0] = xor128_s[0];  s1[1] = xor128_s[1];
  for(int nj=0; nj < 3; nj++) xor128_jump();
  s2[0] = xor128_s[0];  s2[1] = xor128_s[1];

  int nran = 100000;
  double sum1 = 0;      double sum11 = 0;
  double sum2 = 0;      double sum22 = 0;
  double sum12 = 0;
  for(int n=0; n < nran; n++) {
    double ran1 = rand();
    sum1 += ran1;
    sum11 += (ran1-0.5)*(ran1-0.5);
    double ran2 = rand();
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

