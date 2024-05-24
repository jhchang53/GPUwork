#ifndef _inc_XS
#define _inc_XS
/*
	XS structure
*/
struct XSfast {
  double awr;	// for neutron
  double *xsabs;
  double *elas;
  int *nxmus;	// xmu points (size nfast_gps)
  double **xmu;	// equi. prob xmu
  int inel_grps, inel_cdfs;
  double *sig_inel;
  double *cdf_inel,*aniso_inel;
  int n2n_grps,n2n_cdfs;
  double *sig_n2n;
  double *cdf_n2n;
  int n3n_grps,n3n_cdfs;
  double *sig_n3n;
  double *cdf_n3n;
  double *sigtot;
};

struct XStherm {
  int jscat,nbound;
  double *siga,*sigs;
  double *P0,*P1;
  double *sigtot;
};	

struct XSgamma {
  double *phab,*compt,*pair;
  double *sigtot;
};

struct XSgprod {
  double sig_2200,sig_gprod;
  int num_gamline;
  double *g_yield;
  double *g_energy;
};

struct XSkerma {
  double *xkerma_n;
  double *xkerma_g;
};

struct  XS {
  XSfast fast;
  XStherm therm;
  XSgamma gamma;
  XSgprod gprod;
  XSkerma kerma;
};

struct XSL {
  int num_fast_gps, num_therm_gps, num_gam_gps;
  int num_iso;
  double *ev_neut,*ev_gam;
  int *id_fast,*id_therm,*id_gam;
  XS *xs;
};
#endif
