/*
	Mater.h
*/
#include <string>
#include <vector>
#include "XSLIB.h"

using namespace std;

class Mater {
public:
  Mater();
  ~Mater();
  int readMat(const char *matfn);
  int select(const char matname[]);
  void setXSL(XSL *xsl);
  void collect(int m);
  int setupMatLib(std::vector<string> matlist);
private:
  int makeIsoList(int m, int *isolist, double *denlist);
  int collectFast(int m,  int niso, int isolist[], double denlist[],
	double *fast_mac, double *fast_cdf);
  int collectTherm(int m,  int niso, int isolist[], double denlist[],
	double *therm_mac, double *therm_cdf);
  int prlev;
  int nmix;
  char **matname;
  int *nelem;
  int *indmat;
  int *num_pri,*num_sec;
  double *dens;
  XSL *xsl;
  int ng_fast,ng_therm;
  int num_iso;
};
