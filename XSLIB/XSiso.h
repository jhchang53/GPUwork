
#include <stdlib.h>
#include "XS.h"

class XSiso {
public:
  XSiso();
  ~XSiso();
  void setPrlev(int prlev);
  void read(FILE *F, int id_fast, int id_therm, int id_gam,
    int nfast_gps, int ntherm_gps, int ngam_gps);
  XS collect();
private:
  int prlev;
  int id_fast,id_therm,id_gam;
  int nfast_gps,ntherm_gps,ngam_gps;
  int ltot_fast,iwa_fast,iwf_fast,iws_fast,lol_fast[5],la_fast[5],ld_fast[5];
  int jscat,nbound;
  int nrr;
  double awr,spot;
  int ltot_gam,iwa_gam,iwf_gam,iws_gam,lol_gam[5],la_gam[5],ld_gam[5];

  double *adum_fast,*adum_gam;
  double *siga,*sigs,*P0,*P1;	// thermal

  int num_gamline;
  double sig_2200,sig_gprod;
  double *g_yield,*g_energy;

  double *xkerma_n,*xkerma_g;

  void readData(FILE *F, int ndata, double *data);
  int microblk(FILE *F, int ngrpx, double *adum,
    int ltot, int iwa, int iwf, int iws, int lol[], int la[], int ld[]);
  void asm_fast(FILE *F, int nfast_gps);
  void asm_thermal(FILE *F, int ntherm_gps);
  void asm_gamma(FILE *F, int ngam_gps);

  XSfast collect_fast();
  XStherm collect_therm();
  XSgamma collect_gamma();
  XSgprod collect_gprod();
  XSkerma collect_kerma();
};

