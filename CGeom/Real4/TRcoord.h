#ifndef _inc_TRcoord
#define _inc_TRcoord

/*
	TRcoord.h - float version
*/
class TRcoord {
public:
  TRcoord();
  ~TRcoord();
  void setPrlev(int prlev);
  void set(const float target[], const float Zb, const float theta, const float phi, const float psi);
  void S2U(const float q[], float p[]);
  void S2Uvec(const float qvec[], float pvec[]);

  void U2S(const float p[], float s[]);
  void U2Svec(const float pvec[], float svec[]);

private:
  int prlev;
  float PI;
  float deg_rad;
  float R[9];
  float Zb;
  float Tx,Ty,Tz;
  float PTx,PTy,PTz;
};
#endif
