/*
	Tally.h
	need atomic operation
*/

class Tally {
public:
  Tally();
  ~Tally();
  int setLattice(const int nxyz[], const double dxyz[], char *uv);
  void clear();
  void TLE(const double p[], const double w[], const char chx, const double dist);
  float *tal;
  void print();
  void print(const float *tval);
private:
  int prlev;
  int nx,ny,nz;
  int nvoxel;
  double dx,dy,dz;
  char *uv;
  double INF,eps;
};
