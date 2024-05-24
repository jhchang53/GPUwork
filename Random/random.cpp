#include <stdio.h>

#ifndef MS_RAND
#define RAND_MAX ((1U << 31) - 1)
#endif

unsigned s1, s2, s3, s4;

float lcg()
{
  s4 = (( 1664525U*s4 +  1013904223U) & RAND_MAX);
  // return s4*2.3283064365e-10;
  return (float)s4 / RAND_MAX;
};


float taus88()
{
  /* Generates numbers between 0 and 1. */
  unsigned long b;
  unsigned c;
  b = (((s1 << 13) ^ s1) >> 19);
  s1 = (((s1 & 4294967294UL) << 12) ^ b);
  b = (((s2 << 2) ^ s2) >> 25);
  s2 = (((s2 & 4294967288UL) << 4) ^ b);
  b = (((s3 << 3) ^ s3) >> 11);
  s3 = (((s3 & 4294967280UL) << 17) ^ b);
  /*  LCGstep  */
  s4 = (1664525U*s4 + 1013904223U); //  & RAND_MAX;
  // c = s4 >> 16;
  return ((float)(s1 ^ s2 ^ s3 ^ s4) * 2.3283064365e-10);
};

double 

int main()
{
  // s1, s2, s3, s4 should > 128
  s1 = 129;
  s2 = 129;
  s3 = 129;
  s4 = 129;
  for(int i=0; i < 100; i++) {
    float x = taus88();
    printf("%x %.4f\n",s4,x);
    // float x =  lcg();
    // printf("%lx %.4f\n",s4,x);
  }
};
