/*
	CGset.h	- float verssion
*/
#include "CG.h"

class CGset {
public:
  CGset();
  ~CGset();
  void setPrlev(int prlev);
  void set(int ncgs, CG **cg);
  float track(const float p[], const float w[]);
  int cgthis();
  int cgnext();
private:
  int prlev;
  int ncgs;
  CG **cg;
  int *isok;
  float *t_min,*t_max;
  float eps;
  int *isin,*isnext;
};
