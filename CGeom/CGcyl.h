/*
	CGcyl.h
*/
#include "CG.h"

class CGcyl : public CG {
public:
  CGcyl();
  ~CGcyl();
  void CYL(const double R, const double H0, const double H1);
  int track(const double p[], const double w[], double trk[]);
private:
  double R,H0,H1;
};
