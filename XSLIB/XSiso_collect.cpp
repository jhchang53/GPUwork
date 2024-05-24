/*
	XSiso_collect.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "XSiso.h"

XS XSiso::collect()
{
  /*  collect data  */
  XS xs;
  xs.fast = collect_fast();
  xs.therm = collect_therm();
  xs.gamma = collect_gamma();
  xs.gprod = collect_gprod();
  xs.kerma = collect_kerma();
  return xs;
};

XSfast XSiso::collect_fast()
{
  if(prlev) printf("=+++id_fast=%d +++=\n",id_fast);
  double *sigtot = new double[nfast_gps];
  double *xsabs = new double[nfast_gps];
  int iad = 0;
  assert(iwa_fast > 0);
  if(iwa_fast > 0) {
    if(prlev) printf("xsab:");
    for(int g=0; g < nfast_gps; g++) {
      if(prlev) printf(" %.5le",adum_fast[iad]);
      xsabs[g] = adum_fast[iad];
      sigtot[g] = xsabs[g];
      iad++;
      if(prlev) {
        if(g%6==5) printf("\n     ");
      }
    }
    if(prlev) printf("\n");
  }
  assert(iwf_fast == 0);
  double *elas = new double[nfast_gps];
  if(iws_fast > 0) {	// elastic xs
    if(prlev) printf("xsel:");
    for(int g=0; g < nfast_gps; g++) {
      if(prlev) printf(" %.5le",adum_fast[iad]);
      elas[g] = adum_fast[iad];
      sigtot[g] += elas[g];
      iad++;
      if(prlev) {
        if(g%6==5) printf("\n     ");
      }
    }
    if(prlev) printf("\n");
  }
  int inel_grps = 0;
  double *sig_inel = nullptr;
  if(la_fast[0] > 0) {	// inelastic xs
    if(prlev > 1) printf("inel:");
    inel_grps = la_fast[0];
    sig_inel = new double[inel_grps];
    for(int g=0; g < la_fast[0]; g++) {
      if(prlev > 1) printf(" %.5le",adum_fast[iad]);
      sig_inel[g] = adum_fast[iad];
      sigtot[g] += sig_inel[g];
      iad++;
      if(prlev > 1) {
       if(g%6==5) printf("\n     ");
      }
    }
    if(prlev > 1) printf("\n");
  }
  int n2n_grps = 0;
  double *sig_n2n = NULL;
  if(la_fast[1] > 0) {  // n2n
    if(prlev) printf(" n2n:");
    n2n_grps = la_fast[1];
    sig_n2n = new double[n2n_grps];
    for(int g=0; g < la_fast[1]; g++) {
      if(prlev) printf(" %.5le",adum_fast[iad]);
      sig_n2n[g] = adum_fast[iad];
      sigtot[g] += sig_n2n[g];
      iad++;
      if(prlev) {
      if(g%6==5) printf("\n     ");
      }
    }
    if(prlev) printf("\n");
  }
  int n3n_grps = 0;
  double *sig_n3n = NULL;
  if(la_fast[4] > 0) {  // n3n
    if(prlev) printf(" n3n:");
    n3n_grps = la_fast[4];
    sig_n3n = new double[n3n_grps];
    for(int g=0; g < la_fast[4]; g++) {
      if(prlev) printf(" %.5le",adum_fast[iad]);
      sig_n3n[g] = adum_fast[iad];
      sigtot[g] += sig_n3n[g];
      iad++;
      if(prlev) {
      if(g%6==5) printf("\n     ");
      }
    }
    if(prlev) printf("\n");
  }
  int inel_cdfs = 0;
  double *cdf_inel = nullptr;
  if(lol_fast[0] > 0) {	// inelastic down cdf
    if(la_fast[0] != inel_grps) {
      printf("*** la_fast[0]=%d inel_grps=%d\n",la_fast[0],inel_grps);
    }
    assert(la_fast[0] == inel_grps);	// we assume cdf is given to all inel
    inel_cdfs = ld_fast[0];
    int lcdf = inel_cdfs+1;
    cdf_inel = new double[inel_grps*lcdf];
    for(int ig=0; ig < la_fast[0]; ig++) {
      for(int jg=0; jg < lcdf; jg++) {
        double cdf = adum_fast[iad];
        cdf_inel[ig*lcdf+jg] = cdf;
        iad++;
      }
    }
  }
  int inel_aniso = 0;
  double *aniso_inel = nullptr;
  if(lol_fast[2] > 0) {      // inelastic anisotropy
    if(prlev) printf(" collect iad=%d at aniso.\n",iad);
    int nsg = ld_fast[0];
    if(prlev) printf("  inel_grps=%d inel_aniso=%d\n",inel_grps,inel_aniso);
    int lcdf = inel_cdfs+1;
    aniso_inel = new double[inel_grps*lcdf];
    assert(la_fast[0] == inel_grps);
    for(int g=0; g < la_fast[2]-1; g++) {
      for(int j=0; j < lcdf; j++) aniso_inel[g*lcdf+j] = 0;
    }
    for(int g=la_fast[2]-1; g < la_fast[0]; g++) {
      for(int j=0; j < lcdf; j++) {
        aniso_inel[g*lcdf+j] = adum_fast[iad];
        iad++;
      }
    }
    if(prlev) printf(" collect iad=%d at end aniso.\n",iad);
  }
  int n2n_cdfs = 0;
  double *cdf_n2n = nullptr;
  if(lol_fast[1] > 0) {      // (n,2n) down scatter cdf  and photon n2n cdf
    if(prlev) printf(" collect iad=%d at n2n cdf\n",iad);
    assert(la_fast[1] == n2n_grps);
    n2n_cdfs = ld_fast[1];
    int ndat = n2n_cdfs+1;
    cdf_n2n = new double[n2n_grps*ndat];
    if(prlev) printf("  ndat=%d\n",ndat);
    for(int ig=0; ig < la_fast[1]; ig++) {
      if(ndat > 0) {
        for(int jg=0; jg < ndat; jg++) {
          cdf_n2n[ig*ndat+jg] = adum_fast[iad];
          iad++;
        }
        if(prlev) {
          printf("n2n(%d):",ig);
          for(int i=0; i < ndat; i++) {
            printf(" %.6le",cdf_n2n[ig*ndat+i]);
            if(i%6==5) printf("\n  ");
          }
          printf("\n");
        }
      }
    }
    if(prlev) printf(" collect iad=%d at end of n2n cdf\n",iad);
  }
  int n3n_cdfs = 0;
  double *cdf_n3n = nullptr;
  if(lol_fast[4] > 0) {      // (n,3n) downscatter cdf
    if(prlev) printf(" collect iad=%d at n3n cdf\n",iad);
    assert(la_fast[4] == n3n_grps);
    n3n_cdfs = ld_fast[4];
    int ndat3 = n3n_cdfs+1;
    cdf_n3n = new double[n3n_grps*ndat3];
    if(prlev) printf(" la[4]=%d ndat3=%d ld[4]=%d\n",la_fast[4],ndat3,ld_fast[4]);

    if(ld_fast[4] > 0) { // to handle exception of 74-W data
      for(int ig=0; ig < la_fast[4]; ig++) {
        for(int jg=0; jg < ndat3; jg++) {
          cdf_n3n[ig*ndat3+jg] = adum_fast[iad];
          iad++;
        }
      }
    }	// special handling
    else {
      for(int ig=0; ig < la_fast[4]; ig++) {
        int ndd = adum_fast[iad]+0.001;  iad++;
        cdf_n3n[ig*ndat3] = 1.0;
      }
    }
    if(prlev) {
    printf(" cdf_n3n n3n_grps=%d ndat3=%d:\n",n3n_grps,ndat3);
    for(int ig=0; ig < n3n_grps; ig++) {
      printf("%d:",ig);
      for(int jg=0; jg < ndat3; jg++) {
        printf(" %.5le",cdf_n3n[ig*ndat3+jg]);
        if(jg%6 == 5) printf("\n  ");
      }
      printf("\n");
    }
    printf(" collect iad=%d at end of n3n cdf\n",iad);
    }
  }
  /*  elastic eq area xmu  */
  /*  it is mandatory  */
  if(prlev) {
  printf("  collect iad=%d at xmu.\n",iad);
  printf("  nfast_gps=%d\n",nfast_gps);
  }
  int *nxmus = new int[nfast_gps];
  double **xmu = new double*[nfast_gps];
  int iadbeg = iad;
  for(int ig=0; ig < nfast_gps; ig++) {
    int ndim = adum_fast[iad]+0.001;  iad++;
    nxmus[ig] = ndim;
    xmu[ig] = nullptr;
    if(ndim > 0) {
      xmu[ig] = new double[ndim];
      for(int k=0; k < ndim; k++) {
        xmu[ig][k] = adum_fast[iad];
        iad++;
      }
    }
  }
  if(prlev) printf(" collect ltot_fast=%d iad=%d\n",ltot_fast,iad);
  assert(ltot_fast == iad);

  /*  make structure  */
  XSfast xsfast;
  xsfast.awr = awr;
  xsfast.xsabs = xsabs;
  xsfast.elas = elas;
  if(prlev) printf("=== id_fast=%d ===\n",id_fast);
  xsfast.nxmus = nxmus;
  xsfast.xmu = xmu;
  /*  inelastic  */
  xsfast.inel_grps = inel_grps;
  xsfast.sig_inel = sig_inel;
  xsfast.inel_cdfs = inel_cdfs;
  xsfast.cdf_inel = cdf_inel;
  xsfast.aniso_inel = aniso_inel;
  /*  n2n  */
  xsfast.n2n_grps = n2n_grps;
  xsfast.sig_n2n = sig_n2n;
  xsfast.n2n_cdfs = n2n_cdfs;
  xsfast.cdf_n2n = cdf_n2n;
  /*  n3n  */
  xsfast.n3n_grps = n3n_grps;
  xsfast.sig_n3n = sig_n3n;
  xsfast.n3n_cdfs = n3n_cdfs;
  xsfast.cdf_n3n = cdf_n3n;
  /*  total  */
  xsfast.sigtot = sigtot;
  return xsfast;
};	// XSiso::collect_fast

XStherm XSiso::collect_therm()
{
  if(prlev) printf("  collect therm grps =%d\n",ntherm_gps);
  XStherm xstherm;
  xstherm.jscat = jscat;
  xstherm.nbound = nbound;

  if(siga == nullptr) {	// when no thermal data present
    xstherm.siga = nullptr;
    xstherm.sigtot = nullptr;
    return xstherm;
  }
  double *sigtot = new double[ntherm_gps];
  for(int g=0; g < ntherm_gps; g++) {
    sigtot[g] = siga[g]+sigs[g];
  }
  xstherm.siga = siga;
  xstherm.sigs = sigs;
  xstherm.P0 = P0;
  xstherm.P1 = P1;
  xstherm.sigtot = sigtot;
  return xstherm;
};	// XSiso::collect_therm

XSgamma XSiso::collect_gamma()
{
  // sig_phab,*sig_comp,*sig_pair
  if(prlev) printf("=+++id_gam=%d +++=\n",id_gam);
  assert(ngam_gps > 0);
  double *sigtot = new double[ngam_gps];
  double *phab = new double[ngam_gps];
  int iad = 0;
  assert(iwa_gam > 0);
  if(prlev) printf("phab:");
  for(int g=0; g < ngam_gps; g++) {
    if(prlev) printf(" %.5le",adum_gam[iad]);
    phab[g] = adum_gam[iad];
    sigtot[g] = phab[g];
    iad++;
    if(prlev) {
      if(g%6==5) printf("\n     ");
    }
  }
  if(prlev) printf("\n");
  assert(iwf_gam == 0);
  assert(iws_gam > 0);
  double *compt = new double[ngam_gps];
  if(prlev) printf("comp:");
  for(int g=0; g < ngam_gps; g++) {
    if(prlev) printf(" %.5le",adum_gam[iad]);
    compt[g] = adum_gam[iad];
    sigtot[g] += compt[g];
    iad++;
    if(prlev) {
    if(g%6==5) printf("\n     ");
    }
  }
  if(prlev) printf("\n");
  assert(la_gam[0] == 0);
  if(prlev) printf("la[1]=%d\n",la_gam[1]);
  /*  we enforce pair production from group 34 */
  double *pair = new double[34];
  if(prlev) printf("pair:");
  for(int g=0; g < la_gam[1]; g++) {
    if(prlev) printf(" %.5le",adum_gam[iad]);
    pair[g] = adum_gam[iad]; iad++;
    sigtot[g] += pair[g];
    if(prlev) {
      if(g%6==5) printf("\n     ");
    }
  }
  if(prlev) printf("\n");
  /*  for la_gam[1] < 34 case. 26-Fe-Nat  */
  for(int g=la_gam[1]; g < 34; g++) pair[g] = 0.0;
  XSgamma gamma;
  gamma.phab = phab;
  gamma.compt = compt;
  gamma.pair = pair;
  gamma.sigtot = sigtot;
  return gamma;
};	// :XSiso::collect_gamma

XSgprod XSiso::collect_gprod()
{
  XSgprod gprod;
  gprod.num_gamline = num_gamline;
  gprod.sig_2200 = sig_2200;
  gprod.sig_gprod = sig_gprod;
  gprod.g_yield = g_yield;
  gprod.g_energy = g_energy;
  return gprod;
};	// XSiso::collect_gprod

XSkerma XSiso::collect_kerma()
{
  XSkerma kerma;
  kerma.xkerma_n = xkerma_n;
  kerma.xkerma_g = xkerma_g;
  return kerma;
};	// XSiso::collect_kerma
