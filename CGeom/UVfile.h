/*
	UVfile.h
*/
class UVfile {
public:
  UVfile();
  ~UVfile();
  void open(const char *uvpath, const int nxyz[], const double dxyz[]);
  char *UV() { return uv; };
private:
  int nx,ny,nz;
  double dx,dy,dz;
  char *uv;
};
