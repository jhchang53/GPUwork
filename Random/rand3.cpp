#include <stdio.h>

#define a 1
#define b 3
#define c 10

unsigned long x,y,z,w;

unsigned long rand()
{
  unsigned long t;
  t = (x^(x<<a)); x=y; y=z; z=w;
  return w=(w^(w>>c)) ^ (t^(t>>b));
};

int main()
{
  x = 97;
  y = 113;
  for(int i=0; i < 100; i++) rand();

  for(int i=0; i < 20; i++) {
    unsigned long r = rand();
    double ran = r*2.3283064365387e-10;
    printf("%.4lf %ld \n",ran,r);
  }
};
