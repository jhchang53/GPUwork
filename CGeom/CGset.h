/*
	CGset.h
*/
#include "CG.h"

class CGset {
public:
  CGset();
  ~CGset();
  void setPrlev(int prlev);
  void set(const int ncgs, CG **cg);
  double track(const double p[], const double w[]);
  int cgthis();
  int cgnext();
private:
  int prlev;
  int ncgs;
  CG **cg;
  int *isok;
  double *t_min,*t_max;
  double eps;
  int *isin,*isnext;
};
