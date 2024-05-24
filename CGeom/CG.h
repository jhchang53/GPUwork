#ifndef _inc_CG
#define _inc_CG
class CG {
public:
  CG();
  ~CG();
  virtual void CYL(const double R, const double H0, const double H1);
  virtual void TRC(const double R0, const double R1, const double Z0, const double Z1);
  virtual void RPP(const double bmin[], const double bmax[]);
  virtual int track(const double p[], const double w[], double trk[]);
protected:
  int prlev;
  double eps,INF;
};
#endif
