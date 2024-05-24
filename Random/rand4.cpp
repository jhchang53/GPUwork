#include <stdio.h>


unsigned long rand()
{
  static unsigned long x=123456789,y=362436069,z=521288629,
	w=88675123,v=5783321, d=6615241;
  unsigned long t;
  t = (x^(x>>2));  x=y; y=z; z=w; w=v;
  v = (v^(v<<4))^(t^(t<<1));
  return (d+=362437)+v;
};

int main()
{
  // xnorm = 1/18446744073709551615  max of unsigned long
  double xnorm = 5.42101086242752217e-20;
  for(int i=0; i < 100; i++) rand();
  int nran = 5000;
  double xsum = 0;
  for(int i=0; i < nran; i++) {
    
    unsigned long r = rand();
    double ran =  r * xnorm;
    printf("%.4lf %lu \n",ran,r);
    xsum += ran;
  }
  printf("xavg=%.5lf\n",xsum/nran);
  printf("%.20le\n",1.0/18446744073709551615.0);
};
