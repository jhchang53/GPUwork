/*	CGtrc.h  */
#include "CG.h"

class CGtrc : public CG {
public:
  CGtrc();
  ~CGtrc();
  void TRC(const double R0, const double R1, const double Z0, const double Z1);
  int track(const double p[], const double w[], double trk[]);
public:
  double R0,R1,Z0,Z1;
  double Vz,R0sq,R1sq;
  double H,gamsq;
};
