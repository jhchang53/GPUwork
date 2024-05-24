#ifndef _inc_XSLIB
#define _inc_XSLIB
/*
	XSLIB.h
*/
#include "XS.h"

class XSLIB {
public:
  XSLIB();
  ~XSLIB();	// should not destruct until XSL is used
  XSL *openlib(const char *libfn);
  void printXS(XS xs);
private:
  int prlev;
  int num_fast_gps, num_therm_gps, num_gam_gps, num_iso;
  int *id_fast,*id_therm,*id_gam;
  int ngp_neut,ngp_gam;
  double *ev_neut,*ev_gam;

  void printXSfast(XSfast xsfast);
  void printXStherm(XStherm xstherm);

};
#endif
