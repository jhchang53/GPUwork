#ifndef _inc_Scatter
#define _inc_Scatter
/*
        Scatter.h
*/
class Scatter {
public:
  Scatter();
  ~Scatter();
  void setPrlev(int prlev);
  double isotropic(const double Amass, const double E, double w[]);
  double P1(const double p0, const double p1, const double w[], double wvec[]);
  double Compton(const double En);
private:
  void rotate(const double xmu, double k[], double wvec[]);

  int prlev;
  double PI;
};
#endif

