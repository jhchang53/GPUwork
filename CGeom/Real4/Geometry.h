/*
	Geometry.h - float version
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
  Ray track(const float p[], const float w[]);
  void setCG(int ncgs, CG **cg);
  Ray CGtrack(float CGp[], float CGw[]);
  void setLattice(int nxyz[], float dxyz[], char *uv);
  Ray UVtrack(float UVp[], float UVw[], char chx);
  int InUV() { return inUV; };
private:
  int prlev;
  int inUV;	// flag to check
  char chx;

  TRcoord *trc;
  float CGp[3],CGw[3];
  float UVp[3],UVw[3];

  int ncgs;
  CG **cg;
  float INF;
  int *isok;	// flag for line intersection
  float *t_min,*t_max;	// track length for line intersection
  
  Box *box;
  Vox *vox;

  int nextzone;
};
