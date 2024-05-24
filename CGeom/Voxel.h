/*
	Voxel.h
*/
class Voxel {
public:
  Voxel();
  ~Voxel();
  void setLat2D(double dx, double x0, int nx, double dy, double y0, int ny, char*  uv);
  int traverse2D(double px, double py, double wx, double wy);
  void setLat3D(double dxyz[], double xyz0[], int nxyz[], char* uv);
  int traverse3D(double pxyz[], double wxyz[]);
private:
  int nx,ny,nz;
  double x0,y0,z0,dx,dy,dz,xb,yb,zb;
  char *uv;
};
