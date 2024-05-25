/*
	Mater.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <vector>
#include "Mater.h"

using namespace std;

Mater::Mater()
{
  prlev = 0;
  nmix = 0;
  matname = nullptr;
  nelem = nullptr;
  indmat = nullptr;
  num_pri = nullptr;
  num_sec = nullptr;
  dens = nullptr;
};

Mater::~Mater()
{
  for(int m=0; m < nmix; m++) delete [] matname[m];
  delete [] matname;
  delete [] nelem;
  delete [] indmat;
  delete [] num_pri;
  delete [] num_sec;
  delete [] dens;
};

int Mater::readMat(const char *matfn)
{
  FILE *F = fopen(matfn,"r");
  if(F == NULL) {
    printf("*** material data file %s not found.\n",matfn);
    exit(0);
  }
  char buf[100];
  /*  skip comments  */
  fgets(buf,100,F);
  while(buf[0] == '#') fgets(buf,100,F);
  nmix = atoi(buf);
  if(prlev) printf(" nmix=%d\n",nmix);
  matname = new char*[nmix];
  nelem = new int[nmix];
  indmat = new int[nmix+1];
  char matnam[80];
  int maxline = 1000;
  num_pri = new int[maxline];	// primary isotope number (1 bias)
  dens = new double[maxline];
  num_sec = new int[maxline];
  int ind = 0;
  for(int m=0; m < nmix; m++) {
    fgets(buf,100,F);
    while(buf[0] == '#') fgets(buf,100,F);
    char *word = strtok(buf," \n");
    strcpy(matnam,word);
    matname[m] = new char[strlen(matnam)+1];
    strcpy(matname[m],matnam);
    fgets(buf,100,F);
    while(buf[0] == '#') fgets(buf,100,F);
    int nele = atoi(buf);
    nelem[m] = nele;
    for(int n=0; n < nele; n++) {
      fgets(buf,100,F);
      while(buf[0] == '#') fgets(buf,100,F);
      char striso[7];
      for(int i=0; i < 5; i++) striso[i] = buf[i];
      striso[5] = '\0';
      num_pri[ind+n] = atoi(striso);
      char strden[16];
      for(int i=0; i < 15; i++) {
        strden[i] = buf[5+i];
      }
      strden[15] = '\0';
      char strd[16];
      int j=0;
      for(int i=0; i < 15; i++) {
        char ch = strden[i];
        if(ch != ' ') {
          strd[j] = ch;
          j++;
        }
      }
      strd[j] = '\0';
      dens[ind+n] = atof(strd);
      char strsca[6];
      for(int i=0; i < 5; i++) strsca[i] = buf[20+i];
      strsca[5] = '\0';
      num_sec[ind+n] = atoi(strsca);
    }
    indmat[m] = ind;
    ind += nele;
  }
  int indlast = ind;
  indmat[nmix] = indlast;
  if(prlev > 1) {
    ind = 0;
    for(int m=0; m < nmix; m++) {
      int nele = nelem[m];
      printf("%d: %d %s\n",m,nelem[m],matname[m]);
      for(int n=0; n < nele; n++) {
        printf("  %d %.5le %d\n",num_pri[ind+n],dens[ind+n],num_sec[ind+n]);
      }
      ind += nele;
    }
  }
  return nmix;
};	// Mater::readMat

int Mater::select(const char name[])
{
  int mok = -1;
  for(int m=0;(mok < 0) && (m < nmix); m++) {
    if(!strcmp(name,matname[m])) {
      mok = m;
      break;
    }
  }
  if(prlev) {
    if(mok < 0) printf(" %s not found.\n",name);
    else printf(" %s is %s\n",name,matname[mok]);
  }
  return mok;
};	// Mater::select

void Mater::setXSL(XSL *xsl_i)
{
  xsl = xsl_i;
  ng_fast  = xsl->num_fast_gps;
  ng_therm = xsl->num_therm_gps;
  num_iso = xsl->num_iso;
};	// 

void Mater::collect(int m)
{
  if(prlev) printf("%d: %d %s\n",m,nelem[m],matname[m]);
  /*  make full isotope list  */
  int nele = nelem[m];
  int *isolist = new int[2*nele];
  double *denlist = new double[2*nele];
  int nisos = makeIsoList(m, isolist,denlist);
  if(prlev) {
    for(int i=0; i < nisos; i++) printf(" %5d  %.6le\n",isolist[i],denlist[i]);
  }
  double *fast_mac = new double[ng_fast];
  double *fast_cdf = new double[ng_fast*nisos];
  int ok_fast_cdf =collectFast(m,nisos,isolist,denlist,  fast_mac,fast_cdf);
  double *therm_mac = new double[ng_therm];
  double *therm_cdf = new double[ng_therm*nisos];
  int ok_therm_cdf = collectTherm(m, nisos,isolist,denlist, therm_mac,therm_cdf);
  delete [] fast_mac;
  delete [] fast_cdf;
  delete [] therm_mac;
  delete [] therm_cdf;

  delete [] isolist;
  delete [] denlist;
};	// Mater::collect

int Mater::makeIsoList(int m, int *isolist, double *denlist)
{
  if(prlev) printf("%d: %d %s\n",m,nelem[m],matname[m]);
  /*  make full isotope list  */
  int nele = nelem[m];
  int ind = indmat[m];
  int j = 0;
  for(int n=0; n < nele; n++) {
    if(prlev) printf("  %d %.5le %d\n",num_pri[ind+n],dens[ind+n],num_sec[ind+n]);
    /*  primary isotope is always present  */
    isolist[j] = num_pri[ind+n] - 1; // one biad
    denlist[j] = dens[ind+n];
    j++;
    int numsec = num_sec[ind+n];
    if(numsec > 0) {
      int nbound = xsl->xs[numsec-1].therm.nbound;
      isolist[j] = numsec-1;
      denlist[j] = dens[ind+n]/nbound;
      j++;
      assert(nbound > 0);
    }
  }
  int nisos = j;
  return nisos;
};      // Mater::makeIsoList


int Mater::collectFast(int m, int niso, int isolist[], double denlist[],
  double *fast_mac, double *fast_cdf)
{
  /*  macroscopic total cross sections for fast group */
  for(int g=0; g < ng_fast; g++) fast_mac[g] = 0;
  for(int n=0; n < niso; n++) {
    int iso = isolist[n];
    double den = denlist[n];
    XSfast fast = xsl->xs[iso].fast;
    for(int g=0; g < ng_fast; g++) fast_mac[g] += den*fast.sigtot[g];
  }
  if(prlev > 1) {
    printf("fast:");
    for(int g=0; g < ng_fast; g++) {
      printf(" %.5le",fast_mac[g]);
      if(g%6 == 5) printf("\n  ");
    }
    printf("\n");
  }
  /*  chek if total fast is zero  */
  double sumtot = 0.0;
  for(int g=0; g < ng_fast; g++) sumtot += fast_mac[g];
  if(sumtot <= 0.0) return 0;
  /*  prepare fast cdf */
  for(int g=0; g < ng_fast; g++) {
    double summac = 0.0;
    for(int n=0; n < niso; n++) {
      int iso = isolist[n];
      double den = denlist[n];
      XSfast fast = xsl->xs[iso].fast;
      summac += den*fast.sigtot[g];
      fast_cdf[g*niso+n] = summac/fast_mac[g];
    }
  }
  if(prlev > 1) {
    printf(" iso:");
    for(int n=0; n < niso; n++) printf(" %8d",isolist[n]);
    printf("\n");
    for(int g=0; g < ng_fast; g++) {
      printf("fastcdf%2d:",g);
      for(int n=0; n < niso; n++) printf(" %.6lf",fast_cdf[g*niso+n]);
      printf("\n");
    }
  }
  return 1;
};	// Mater::collectFast

int Mater::collectTherm(int m, int niso, int isolist[], double denlist[],
  double *therm_mac, double *therm_cdf)
{
  /*  macroscopic total cross sections for therm group */
  for(int g=0; g < ng_therm; g++) therm_mac[g] = 0;
  for(int n=0; n < niso; n++) {
    int iso = isolist[n];
    double den = denlist[n];
    XStherm therm = xsl->xs[iso].therm;
    /*  thermal cross section can be missing  */
    if(therm.sigtot != nullptr) 
    for(int g=0; g < ng_therm; g++) therm_mac[g] += den*therm.sigtot[g];
  }
  if(prlev > 1) {
    printf("therm:");
    for(int g=0; g < ng_therm; g++) {
      printf(" %.5le",therm_mac[g]);
      if(g%6 == 5) printf("\n  ");
    }
    printf("\n");
  }
  /*  chek if total thermal is zero  */
  double sumtot = 0.0;
  for(int g=0; g < ng_therm; g++) sumtot += therm_mac[g];
  if(sumtot <= 0.0) return 0;

  /*  prepare therm cdf */
  for(int g=0; g < ng_therm; g++) {
    double summac = 0.0;
    for(int n=0; n < niso; n++) {
      int iso = isolist[n];
      double den = denlist[n];
      XStherm therm = xsl->xs[iso].therm;
      if(therm.sigtot != nullptr) summac += den*therm.sigtot[g];
      therm_cdf[g*niso+n] = summac/therm_mac[g];
    }
  }
  if(prlev > 1) {
    printf(" iso:");
    for(int n=0; n < niso; n++) printf(" %8d",isolist[n]);
    printf("\n");
    for(int g=0; g < ng_therm; g++) {
      printf("thermcdf%2d:",g);
      for(int n=0; n < niso; n++) printf(" %.6lf",therm_cdf[g*niso+n]);
      printf("\n");
    }
  }
  return 1;
};	//  Mater::collectTherm

int Mater::setupMatLib(std::vector<string> matlist)
{
  nmat = matlist.size();
  ind_iso = new int[nmat+1];
  iso_vec.clear();
  ind_fmac = new int[nmat+1];
  ind_fcdf = new int[nmat+1];
  fast_mac_vec.clear();
  fast_cdf_vec.clear();
  ind_tmac = new int[nmat+1];
  ind_tcdf = new int[nmat+1];
  therm_mac_vec.clear();
  therm_cdf_vec.clear();

  int indiso = 0;
  int indfm = 0;
  int indfc = 0;
  int indtm = 0;
  int indtc = 0;
  for(int mat=0; mat < nmat; mat++) {
    std::string name = matlist[mat];
    if(prlev) printf(" [%s]\n",name.c_str());
    int m = select(name.c_str());
    /*  make full isotope list  */
    int nele = nelem[m];
    if(prlev) printf("%d: %d %s\n",m,nelem[m],matname[m]);
    int *isolist = new int[2*nele];
    double *denlist = new double[2*nele];
    int nisos = makeIsoList(m, isolist,denlist);
    if(prlev) {
      for(int i=0; i < nisos; i++) printf(" %5d  %.6le\n",isolist[i],denlist[i]);
    }
    ind_iso[mat] = indiso;
    std::vector<int> isvec(isolist,isolist+nisos);
    iso_vec.insert(std::end(iso_vec), std::begin(isvec),std::end(isvec));

    indiso += nisos;
    double *fmac = new double[ng_fast];
    double *fcdf = new double[ng_fast*nisos];
    int ok_fast_cdf =collectFast(m,nisos,isolist,denlist,  fmac,fcdf);
    ind_fmac[mat] = indfm;
    ind_fcdf[mat] = indfc;
    if(ok_fast_cdf) {
      std::vector<double> fmac_vec(fmac,fmac+ng_fast);
      fast_mac_vec.insert(std::end(fast_mac_vec),
        std::begin(fmac_vec),std::end(fmac_vec));
      indfm += ng_fast;
      std::vector<double> fcdf_vec(fcdf,fcdf+ng_fast*nisos);
      fast_cdf_vec.insert(std::end(fast_cdf_vec),
        std::begin(fcdf_vec),std::end(fcdf_vec));
      indfc += ng_fast*nisos;
    }
    double *tmac = new double[ng_therm];
    double *tcdf = new double[ng_therm*nisos];
    int ok_therm_cdf = collectTherm(m, nisos,isolist,denlist, tmac,tcdf);
    ind_tmac[mat] = indtm;
    ind_tcdf[mat] = indtc;
    if(ok_therm_cdf) {
      std::vector<double> tmac_vec(tmac,tmac+ng_therm);
      therm_mac_vec.insert(std::end(therm_mac_vec),
        std::begin(tmac_vec),std::end(tmac_vec));
      ind_tmac[mat] = indtm;
      indtm += ng_therm;
      std::vector<double> tcdf_vec(tcdf,tcdf+ng_therm*nisos);
      therm_cdf_vec.insert(std::end(therm_cdf_vec),
        std::begin(tcdf_vec),std::end(tcdf_vec));
      indtc += ng_therm*nisos;
    }
    delete [] fmac;
    delete [] fcdf;
    delete [] tmac;
    delete [] tcdf;

    delete [] isolist;
    delete [] denlist;
  }	// for mat
  ind_iso[nmat] = indiso;
  ind_fmac[nmat] = indfm;
  ind_fcdf[nmat] = indfc;
  ind_tmac[nmat] = indtm;
  ind_tcdf[nmat] = indtc;

  assert(indfm == fast_mac_vec.size());
  return 0;
};	// Mater::setupMatLib

void Mater::print(int look)
{
  printf("fast:");
  for(int i=ind_fmac[look]; i < ind_fmac[look+1]; i++) {
    printf(" %.5le",fast_mac_vec[i]);
    if(i%6==5) printf("\n  ");
  }
  printf("\n");
  int numiso = ind_iso[look+1]-ind_iso[look];
  printf("fast_cdf:");
  for(int i=ind_iso[look]; i < ind_iso[look+1]; i++) {
    printf("%5d       ",iso_vec[i]);
  }
  printf("\n");
  printf("  ");
  int j = 0;
  for(int i=ind_fcdf[look]; i < ind_fcdf[look+1]; i++) {
    printf(" %.5le",fast_cdf_vec[i]);
    if(j%numiso == numiso-1) printf("\n  ");
    j++;
  }
  printf("\n");
  printf("therm:");
  j = 0;
  for(int i=ind_tmac[look]; i < ind_tmac[look+1]; i++) {
    printf(" %.5le",therm_mac_vec[i]);
    if(j%6==5) printf("\n  ");
    j++;
  }
  printf("\n");
  printf("therm cdf:");
  for(int i=ind_iso[look]; i < ind_iso[look+1]; i++) {
    printf("%4d        ",iso_vec[i]);
  }
  printf("\n");
  j = 0;
  printf("  ");
  for(int i=ind_tcdf[look]; i < ind_tcdf[look+1]; i++) {
    printf(" %.5le",therm_cdf_vec[i]);
    if(j%numiso == numiso-1) printf("\n  ");
    j++;
  }
  printf("\n");

};	// Mater::print
