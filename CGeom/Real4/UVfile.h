/*
	UVfile.h - float version
*/
class UVfile {
public:
  UVfile();
  ~UVfile();
  void open(const char *uvpath, const int nxyz[], const float dxyz[]);
  char *UV() { return uv; };
private:
  int nx,ny,nz;
  float dx,dy,dz;
  char *uv;
};
