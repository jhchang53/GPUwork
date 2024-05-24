/*
	Box.h
	ref) D.Eberly, Intgersection of a Line and a Box, Geometric Tools, 2018
*/
class Box {
public:
  Box();
  ~Box();
  void set(double bmin[], double bmax[]);
  int isInside(double P[]);
  int DoRayQuery(const double P[], const double D[], double t[2]);

private:
  int Clip(double denom, double numer, double &t0, double &t1);
  int DoLineQuery(const double P[], const double D[], double t[2]);

  int prlev;
  double xmin[3],xmax[3];
  double box_C[3],box_e[3];
  double INF;
};
