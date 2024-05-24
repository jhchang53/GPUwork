#include <stdio.h>

unsigned long x,y;

unsigned long rand()
{
  unsigned long t;
  t = (x^(x<<2)); x=y; return y=(y^(y>>4)) ^ (t^(t>>1));
};

int main()
{
  x = 97;
  y = 113;
  for(int i=0; i < 20; i++) {
    unsigned long r = rand();
    printf("%ld \n",r);
  }
};
