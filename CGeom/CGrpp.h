/*
	CGrpp.h
	ref) D.Eberly, Intgersection of a Line and a Box, Geometric Tools, 2018
*/
#include "CG.h"

class CGrpp : public CG {
public:
  CGrpp();
  ~CGrpp();
  void RPP(const double bmin[], const double bmax[]);
  int track(const double p[], const double w[], double trk[]);
private:
  int Clip(double denom, double numer, double &t0, double &t1);

  int prlev;
  double box_C[3],box_e[3];
};
