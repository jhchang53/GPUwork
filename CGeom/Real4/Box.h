/*
	Box.h	- float version
	ref) D.Eberly, Intgersection of a Line and a Box, Geometric Tools, 2018
*/
class Box {
public:
  Box();
  ~Box();
  void set(const float bmin[], const float bmax[]);
  int isInside(const float P[]);
  int DoRayQuery(const float P[], const float D[], float t[2]);

private:
  int Clip(float denom, float numer, float &t0, float &t1);
  int DoLineQuery(const float P[], const float D[], float t[2]);

  int prlev;
  float xmin[3],xmax[3];
  float box_C[3],box_e[3];
  float INF;
};
