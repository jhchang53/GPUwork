/*
	Tally.h	- float version
	need atomic operation
*/

class Tally {
public:
  Tally();
  ~Tally();
  int setLattice(const int nxyz[], const float dxyz[], char *uv);
  void clear();
  void TLE(const float p[], const float w[], char chx, float dist);
  float *tal;
  void print();
  void print(float *tval);
private:
  int prlev;
  int nx,ny,nz;
  int nvoxel;
  float dx,dy,dz;
  char *uv;
  float INF,eps;
};
