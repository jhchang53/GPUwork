/*
	Geometry.h
*/
#include "CG.h"
#include "Box.h"
#include "Vox.h"
#include "Ray.h"
#include "TRcoord.h"

class Geometry {
public:
  Geometry();
  ~Geometry();
  void setPrlev(int prlev);
  void setTR(TRcoord *trc);
  Ray track(const double p[], const double w[]);
  void setCG(int ncgs, CG **cg);
  Ray CGtrack(double CGp[], double CGw[]);
  void setLattice(int nxyz[], double dxyz[], char *uv);
  Ray UVtrack(double UVp[], double UVw[], char chx);
  int InUV() { return inUV; };
private:
  int prlev;
  int inUV;	// flag to check
  char chx;

  TRcoord *trc;
  double CGp[3],CGw[3];
  double UVp[3],UVw[3];

  int ncgs;
  CG **cg;
  double INF;
  int *isok;	// flag for line intersection
  double *t_min,*t_max;	// track length for line intersection
  
  Box *box;
  Vox *vox;

  int nextzone;
};
