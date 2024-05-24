/*
	CGrpp.h
	ref) D.Eberly, Intgersection of a Line and a Box, Geometric Tools, 2018
*/
#include "CG.h"

class CGrpp : public CG {
public:
  CGrpp();
  ~CGrpp();
  void RPP(const float bmin[], const float bmax[]);
  int track(const float p[], const float w[], float trk[]);
private:
  int Clip(float denom, float numer, float &t0, float &t1);

  int prlev;
  float box_C[3],box_e[3];
};
