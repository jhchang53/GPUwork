#include <stdio.h>

#define a 7
#define b 13
#define c 6

unsigned long x,y,z;

unsigned long rand()
{
  unsigned long t;
  t = (x^(x<<a)); x=y; y=z;
  return z=(z^(z>>c)) ^ (t^(t>>b));
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
