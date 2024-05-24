/*
	Vox.h - float version
*/
#include "Ray.h"

class Vox {
public:
  Vox();
  ~Vox();
  void setLattice(const int nxyz[], const float dxyz[], char* uv);
  char getChar(const float pxyz[]);
  Ray track(const float pxyz[], const float wxyz[], char chx);
private:
  int prlev;
  int nx,ny,nz;
  float INF,eps;
  float dx,dy,dz;
  char *uv;
};
