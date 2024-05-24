#ifndef _inc_TRcoord
#define _inc_TRcoord

/*
	TRcoord.h
*/
class TRcoord {
public:
  TRcoord();
  ~TRcoord();
  void setPrlev(int prlev);
  void set(const double target[], const double Zb, const double theta, const double phi, const double psi);
  void S2U(const double q[], double p[]);
  void S2Uvec(const double qvec[], double pvec[]);

  void U2S(const double p[], double s[]);
  void U2Svec(const double pvec[], double svec[]);

private:
  int prlev;
  double PI;
  double deg_rad;
  double R[9];
  double Zb;
  double Tx,Ty,Tz;
  double PTx,PTy,PTz;
};
#endif
