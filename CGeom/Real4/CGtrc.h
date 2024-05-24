/*	CGtrc.h  */
#include "CG.h"

class CGtrc : public CG {
public:
  CGtrc();
  ~CGtrc();
  void TRC(const float R0, const float R1, const float Z0, const float Z1);
  int track(const float p[], const float w[], float trk[]);
public:
  float R0,R1,Z0,Z1;
  float Vz,R0sq,R1sq;
  float H,gamsq;
};
