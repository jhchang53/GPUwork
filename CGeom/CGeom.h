/*
	CGeom.h
*/
class CGeom {
public:
  CGeom();
  ~CGeom();

  double trackMax();
  double trackMin();
  int RPP(double xmin[], double xmax[], double p[], double w[]);
  int CYL(double R, double H0, double H1, double p[], double w[]);

  void setLat2D(double dx, double x0, int nx, double y0, int ny, char *uv);
  int locate2D(double px, double py, double wx, double wy);
private:
  int next2D(int ch, int face, double wx, double wy, int ixi, int jyi);
  int prlev;
  double epsilon;
  int nx,ny;
  char *uv;
  double dx,dy,x0,y0;
  int ix,jy;
  double ux,uy;
  double xb,yb;

  double tmin,tmax;
};
