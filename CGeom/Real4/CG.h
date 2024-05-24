#ifndef _inc_CG
#define _inc_CG
class CG {
public:
  CG();
  ~CG();
  virtual void CYL(const float R, const float H0, const float H1);
  virtual void TRC(const float R0, const float R1, const float Z0, const float Z1);
  virtual void RPP(const float bmin[], const float bmax[]);
  virtual int track(const float p[], const float w[], float trk[]);
protected:
  int prlev;
  float eps,INF;
};
#endif
