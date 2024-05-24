/*
	CGrpp.h
*/
#include "CG.h"

class CGrpp : public CG {
public:
  CGrpp();
  ~CGrpp();
  void RPP(double x0, double x1, double y0, double y1, double z0, double z1);
  int track(double p[], double w[], double trk[]);
private:
  double xmin[3],xmax[3];
};
