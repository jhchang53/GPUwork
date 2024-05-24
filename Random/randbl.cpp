#include <stdio.h>

/* ref. L'Ecuyer, P., Math. Comp. 68(224_ 261-269. (1999)  */
unsigned long z1, z2, z3, z4;

double lfsr113()
{
  /* Generates numbers between 0 and 1. */
  unsigned long b;
  b  = (((z1 << 6) ^ z1) >> 13);
  z1 = (((z1 & 4294967294) << 18) ^ b);
  b  = (((z2 << 2) ^ z2) >> 27);
  z2 = (((z2 & 4294967288) << 2) ^ b);
  b  = (((z3 << 13) ^ z3) >> 21);
  z3 = (((z3 & 4294967280) << 7) ^ b);
  b  = (((z4 <<  3) ^ z4) >> 12);
  z4 = (((z4 & 4294967168) << 13) ^ b);
  return ((z1 ^ z2 ^ z3 ^ z4) * 2.3283064365387e-10);
};


int main()
{
  // s1, s2, s3, s4 should > 128
  z1 = 129;
  z2 = 130;
  z3 = 131;
  z4 = 132;
  for(int i=0; i < 100; i++) {
    double x = lfsr113();
    printf("%.4lf\n",x);
  }
};
