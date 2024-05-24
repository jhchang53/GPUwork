/*
	Vox.h
*/
#include "Ray.h"

class Vox {
public:
  Vox();
  ~Vox();
  void setLattice(int nxyz[], double dxyz[], char* uv);
  char getChar(const double pxyz[]);
  Ray track(const double pxyz[], const double wxyz[], char chx);
private:
  int prlev;
  int nx,ny,nz;
  double INF,eps;
  double dx,dy,dz;
  char *uv;
};
