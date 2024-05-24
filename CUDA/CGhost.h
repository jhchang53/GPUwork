/*
	CGhost.h
*/
#define CGBAR 8
int CGrpp_host(double *cgparm, const double xmin[], const double xmax[]);
int CGcyl_host(double *cgparm, const double R, const double H0, const double H1);
int CGtrc_host(double *cgparm, const double R0, const double R1,
  const double Z0, const double Z1);

