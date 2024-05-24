/*
	isotro.cpp
*/
#include <stdio.h>
#include <math.h>
#include "Xor128.h"

int main()
{
  Xor128 *ran = new Xor128();
  for(int nj=0; nj < 3; nj++) ran->jump();

  int nran = 1000;
  for(int n=0; n < nran; n++) {
    /* generate isotropic direction vector.	*/
    /* ref) G.Marsaglia, Choosing a point from the surface of a sphere,   */
    /*		Ann. Math.Stat. 43(2) 645-646, (1972)  */
    double v1 = 2*ran->ranf()-1.0;
    double v2 = 2*ran->ranf()-1.0;
    double S = v1*v1+v2*v2;
    if(S < 1.0) {
      double rtS = sqrt(1.0-S);
      double wx = 2*v1*rtS;
      double wy = 2*v2*rtS;
      double wz = 1.0-2*S;
      printf("%.4lf %.4lf %.4lf   %.4lf\n",wx,wy,wz,wx*wx+wy*wy+wz*wz);
    }
  }
};
