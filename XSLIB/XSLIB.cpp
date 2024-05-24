/*
	XSLIB.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "XSLIB.h"
#include "XSiso.h"

XSLIB::XSLIB()
{
  prlev = 0;
  id_fast = nullptr;
  id_therm = nullptr;
  id_gam = nullptr;
  ev_neut = nullptr;
  ev_gam = nullptr;
};

XSLIB::~XSLIB()
{
  delete [] id_fast;
  delete [] id_therm;
  delete [] id_gam;
  delete [] ev_neut;
  delete [] ev_gam;
};

XSL *XSLIB::openlib(const char *libfn)
{
  XSL *xsl = new XSL();
  FILE *F = fopen(libfn,"r");
  if(F == NULL) {
    printf("*** XS library file %s not found.\n",libfn);
    exit(0);
  }
  /*  library header  */
  char buf[90];
  fgets(buf,90,F);
  if(prlev) printf("%s",buf);
  /*    num_fast_gps num_therm_gps num_gam_gps num_iso  */
  fgets(buf,90,F);
  char *word = strtok(buf," "); num_fast_gps = atoi(word);
  word = strtok(NULL," ");	num_therm_gps = atoi(word);
  word = strtok(NULL," ");      num_gam_gps = atoi(word);
  word = strtok(NULL," ");      num_iso = atoi(word);
  if(prlev) printf("  num_fast_gps=%d num_therm_gps=%d num_gam_gps=%d num_iso=%d\n",
	 num_fast_gps,num_therm_gps,num_gam_gps,num_iso);
  xsl->num_fast_gps = num_fast_gps;
  xsl->num_therm_gps = num_therm_gps;
  xsl->num_gam_gps = num_gam_gps;
  xsl->num_iso = num_iso;
  /*  headers  */
  for(int j=0; j < 5; j++) {
    fgets(buf,90,F);
    if(prlev) printf(" %s",buf);
  }
  /*  ids  */
  id_fast = new int[num_iso];
  id_therm = new int[num_iso];
  id_gam = new int[num_iso];

  for(int iso=0; iso < num_iso; iso++) {
    fgets(buf,90,F);
    word = strtok(buf," ");	id_fast[iso] = atoi(word);
    word = strtok(NULL," ");	id_therm[iso] = atoi(word);
    word = strtok(NULL," ");    id_gam[iso] = atoi(word);
  }
  xsl->id_fast = id_fast;
  xsl->id_therm = id_therm;
  xsl->id_gam = id_gam;

  /*  neutron energy group boundary  */
  ngp_neut = num_fast_gps + num_therm_gps;
  ev_neut = new double[ngp_neut+1];
  int ncards = (ngp_neut+5)/6;
  int g = 0;
  for(int nc=0; nc < ncards; nc++) {
    fgets(buf,90,F);
    for(int j=0; j < 6; j++) {
      if(j==0) word = strtok(buf," ");
      else word = strtok(NULL," ");
      ev_neut[g] = atof(word);
      g++;
      if(g > ngp_neut) break;
    }
  }
  if(prlev) {
    printf("Evn:");
    for(int g=0; g < ngp_neut+1; g++) {
      printf(" %.5le",ev_neut[g]);
      if(g%6 == 5) printf("\n");
    }
    printf("\n");
  }
  xsl->ev_neut = ev_neut;
  /*  gamma energy group boundary  */
  ngp_gam = num_gam_gps;
  ev_gam = new double[ngp_gam+1];
  ncards = (ngp_gam+5)/6;
  g = 0;
  for(int nc=0; nc < ncards; nc++) {
    fgets(buf,90,F);
    for(int j=0; j < 6; j++) {
      if(j==0) word = strtok(buf," ");
      else word = strtok(NULL," ");
      ev_gam[g] = atof(word);
      g++;
      if(g > ngp_gam) break;
    }
  }
  xsl->ev_gam = ev_gam;
  if(prlev) {
    printf("Evg:");
    for(int g=0; g < ngp_gam+1; g++) {
      printf(" %.5le",ev_gam[g]);
      if(g%6 == 5) printf("\n");
    }
    printf("\n");
  }
  XSiso **xsiso = new XSiso*[num_iso];
  int look = -2; // 29; // 19; // 9; // num_iso-4;
  for(int iso=0; iso < num_iso; iso++) {
    xsiso[iso] = new XSiso();
    if(iso == look) xsiso[iso]->setPrlev(1);
    else xsiso[iso]->setPrlev(0);
    xsiso[iso]->read(F, id_fast[iso], id_therm[iso], id_gam[iso], num_fast_gps,
      num_therm_gps, num_gam_gps);
    if(iso == look) {
      if(prlev) printf("+++ id_fast=%d \n",id_fast[iso]);
      XS xs = xsiso[iso]->collect();
      printXS(xs);
    }
  }
  fclose(F);
  if(prlev) printf("=== XSLIB processed. ===\n");
  /*  move to XSL structure  */
  xsl->xs = new XS[num_iso];
  for(int iso=0; iso < num_iso; iso++) {
    xsl->xs[iso] = xsiso[iso]->collect();
    delete xsiso[iso];
  }
  delete [] xsiso;
  return xsl;
};	// XSLIB::openlib

void XSLIB::printXS(XS xs)
{
  // printXSfast(xs.fast);
  printXStherm(xs.therm);
};

void XSLIB::printXSfast(XSfast fast)
{
  printf("= xsabs:");
  for(int g=0; g < num_fast_gps; g++) {
    printf(" %.5le", fast.xsabs[g]);
    if(g%6==5) printf("\n   ");
  }
  printf("\n");
  printf("= elas:");
  for(int g=0; g < num_fast_gps; g++) {
    printf(" %.5le", fast.elas[g]);
    if(g%6==5) printf("\n   ");
  }
  printf("\n");

  if(prlev) {
  printf("= xmu: ngp=%d\n",num_fast_gps);
  for(int g=0; g < num_fast_gps; g++) {
    printf("%d %d:",g,fast.nxmus[g]);
    for(int k=0; k < fast.nxmus[g]; k++) printf(" %.6le",fast.xmu[g][k]);
    printf("\n");
  }
  }
  printf("= inel:");
  for(int g=0; g < fast.inel_grps; g++) {
    printf(" %.5le", fast.sig_inel[g]);
    if(g%6==5) printf("\n   ");
  }
  printf("\n");
  if(prlev) {
     printf("= inel cdf %d:\n",fast.inel_cdfs);
    int lcdf = fast.inel_cdfs + 1;
    for(int ig=0; ig < fast.inel_grps; ig++) {
      printf("%d:",ig);
      for(int jg=0; jg < lcdf; jg++) printf(" %.6le",
	fast.cdf_inel[ig*lcdf+jg]);
      printf("\n");
    }
  }
  if(prlev) {
  printf("= inel aniso %d:\n",fast.inel_cdfs+1);
  int liso = fast.inel_cdfs + 1;
  for(int ig=0; ig < fast.inel_grps; ig++) {
    printf("%d:",ig);
    for(int jg=0; jg < liso; jg++) printf(" %.6le",
        fast.aniso_inel[ig*liso+jg]);
    printf("\n");
  }
  }
  printf("= n2n:");
  for(int g=0; g < fast.n2n_grps; g++) {
    printf(" %.5le", fast.sig_n2n[g]);
    if(g%6==5) printf("\n   ");
  }
  printf("\n");
  if(prlev) {
    printf("= n2n cdf %d:\n",fast.n2n_cdfs);
    int lcdf2 = fast.n2n_cdfs + 1;
    for(int ig=0; ig < fast.n2n_grps; ig++) {
      printf("%d:",ig);
      for(int jg=0; jg < lcdf2; jg++) {
        printf(" %.6le", fast.cdf_n2n[ig*lcdf2+jg]);
        if(jg%6==5) printf("\n   ");
      }
      printf("\n");
    }
  }
  printf("= n3n:");
  for(int g=0; g < fast.n3n_grps; g++) {
    printf(" %.5le", fast.sig_n3n[g]);
    if(g%6==5) printf("\n   ");
  }
  printf("\n");
  printf("= n3n cdf %d:\n",fast.n3n_cdfs);
  int lcdf2 = fast.n3n_cdfs + 1;
  for(int ig=0; ig < fast.n3n_grps; ig++) {
    printf("%d:",ig);
    for(int jg=0; jg < lcdf2; jg++) {
      printf(" %.6le", fast.cdf_n3n[ig*lcdf2+jg]);
      if(jg%6==5) printf("\n   ");
    }
    printf("\n");
  }

};	//  XSLIB::printXSfast

void XSLIB::printXStherm(XStherm therm)
{
  printf(" jscat=%d nbound=%d\n",therm.jscat,therm.nbound);
  if(therm.siga == nullptr) return;
  printf("= siga:");
  for(int g=0; g < num_therm_gps; g++) {
    printf(" %.5le",therm.siga[g]);
    if(g%6 == 5) printf("\n  ");
  }
  printf("\n");
  printf("= sigs:");
  for(int g=0; g < num_therm_gps; g++) {
    printf(" %.5le",therm.sigs[g]);
    if(g%6 == 5) printf("\n  ");
  }
  printf("\n");
  for(int g1=0; g1 < num_therm_gps; g1++) {
    printf("P0 %d:",g1);
    for(int g2=0; g2 < num_therm_gps; g2++) {
      printf(" %.5le",therm.P0[g1*num_therm_gps+g2]);
      if(g2%6 == 5) printf("\n  ");
    }
    printf("\n");
  }
  for(int g1=0; g1 < num_therm_gps; g1++) {
    printf("P1 %d:",g1);
    for(int g2=0; g2 < num_therm_gps; g2++) {
      printf(" %.5le",therm.P1[g1*num_therm_gps+g2]);
      if(g2%6 == 5) printf("\n  ");
    }
    printf("\n");
  }

};	// XSLIB::printXStherm

