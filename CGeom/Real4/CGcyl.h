/*
	CGcyl.h	- floart version
*/
#include "CG.h"

class CGcyl : public CG {
public:
  CGcyl();
  ~CGcyl();
  void CYL(const float R, const float H0, const float H1);
  int track(const float p[], const float w[], float trk[]);
private:
  float R,H0,H1;
};
