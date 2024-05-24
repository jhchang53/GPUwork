/*
	XSiso.h
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "XSiso.h"

XSiso::XSiso()
{
  prlev = 0;
  adum_fast = nullptr;
  adum_gam = nullptr;
  siga = nullptr;
  sigs = nullptr;
  P0 = nullptr;
  P1 = nullptr;
  g_yield = nullptr;
  g_energy = nullptr;
  xkerma_n = nullptr;
  xkerma_g = nullptr;
};

XSiso::~XSiso()
{
  delete [] adum_fast;
  delete [] adum_gam;
  // delete [] g_yield;
  // delete [] g_energy;
};

void XSiso::setPrlev(int prl)
{
  prlev = prl;
};	// XSiso::setPrlev

void XSiso::read(FILE *F, int id_fast_i, int id_therm_i, int id_gam_i,
  int nfast_gps_i, int ntherm_gps_i, int ngam_gps_i)
{
  id_fast = id_fast_i;
  id_therm = id_therm_i;
  id_gam = id_gam_i;

  nfast_gps  = nfast_gps_i;
  ntherm_gps = ntherm_gps_i;
  ngam_gps   = ngam_gps_i;
  char buf[90];
  char *word;
  if(id_fast > 0) {
    if(prlev) printf("== id_fast=%d ==\n",id_fast);
    fgets(buf,90,F);
    if(prlev) printf(" fast = %s",buf);
    asm_fast(F,nfast_gps);
  }
  jscat = 0;
  nbound = 0;
  char strbuf[7];
  if(id_therm != 0) {
    fgets(buf,90,F);
    strncpy(strbuf,buf,6);
    strbuf[6] = '\0';
    jscat = atoi(strbuf);
    strncpy(strbuf,buf+6,6);
    strbuf[6] = '\0';
    nbound = atoi(strbuf);
    asm_thermal(F,ntherm_gps);
  }
  if(id_gam != 0) {
    asm_gamma(F,ngam_gps);
  }
  /*  capture gamma  */
  fgets(buf,90,F);
  if(prlev) printf("   %s",buf);
  fgets(buf,90,F);
  word = strtok(buf," ");	num_gamline = atoi(word);
  word = strtok(NULL," ");	sig_2200 = atof(word);
  word = strtok(NULL," ");      sig_gprod = atof(word);
  if(prlev) printf("  num_gamline=%d sig_2200=%.4le sig_gprod=%.4le\n",
	num_gamline,sig_2200,sig_gprod);
  /*  xkerma_n  */
  fgets(buf,90,F);
  if(prlev) printf("   %s",buf);
  xkerma_n = new double[nfast_gps+ntherm_gps];
  if(prlev) printf(" ngrps=%d\n",nfast_gps+ntherm_gps);
  readData(F,nfast_gps+ntherm_gps,xkerma_n);
  /*  xkerma_g  */
  fgets(buf,90,F);
  if(prlev) printf("   %s",buf);
  xkerma_g = new double[ngam_gps];
  readData(F,ngam_gps,xkerma_g);
  /* gamma yield and energy  */
  g_yield = nullptr;
  g_energy = nullptr;
  if(num_gamline > 0) {
    fgets(buf,90,F);
    if(prlev) printf("   %s",buf);
    g_yield = new double[num_gamline];
    readData(F,num_gamline, g_yield);
    g_energy = new double[num_gamline];
    readData(F,num_gamline, g_energy);
  }
  if(prlev) {
    printf(" === XSiso::read done ===\n");
    if(siga) printf(" siga[0]=%.5le\n",siga[0]);
  }
};	// XSiso::read

void XSiso::readData(FILE *F, int ndata, double *data)
{
  /*  read array of ndata in F12.5 format  */
  char buf[80];
  char word[13];
  int ncard = (ndata+5)/6;
  int n = 0;
  for(int nc=0; nc < ncard; nc++) {
    fgets(buf,80,F);
    for(int j=0; j < 6; j++) {
      for(int i=0; i < 12; i++) word[i] = buf[12*j+i];
      word[12] = '\0';
      data[n] = atof(word);
      n++;
      if(n >= ndata) break;
    }
  }
  assert(n==ndata);
};	// XSiso::readData


int XSiso:: microblk(FILE *F, int ngrpx, double *adum,
  int ltot, int iwa, int iwf, int iws, int lol[], int la[], int ld[])
{
  int iad = 0;
  char buf[90];
  char *word;
  int ncards = (ngrpx+5)/6;
  if(prlev > 1) printf("  ngrpx=%d ncards=%d\n",ngrpx,ncards);
  assert(iwa >= 0);
  if(iwa > 0) {	// neutron absorption xs
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    readData(F, ngrpx, adum+iad);
    iad += ngrpx;
  }
  assert(iwf >= 0);
  if(iwf > 0) {	// neutron fission xs
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    readData(F, ngrpx, adum+iad);
    iad += ngrpx;
  }
  assert(iws >= 0);
  if(iws > 0) { // elastic xs
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    readData(F, ngrpx, adum+iad);
    iad += ngrpx;
  }
  if(la[0] > 0) {	// inelastic xs
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    readData(F, la[0], adum+iad);
    iad += la[0];
  }
  if(la[1] > 0) {       // n2n xs
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    readData(F, la[1], adum+iad);
    iad += la[1];
  }
  if(la[4] > 0) {       // n3n xs
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    readData(F, la[4], adum+iad);
    iad += la[4];
  }
  /*  down scatter probability  */
  if(lol[0] > 0) {	// inelastic down cdf
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    double cdf[100];	// cdf size is 1 large covering from [0.0, 1.0]
    for(int ig=0; ig < la[0]; ig++) {
      fgets(buf,90,F);
      word = strtok(buf," ");   int igr = atoi(word);
      word = strtok(NULL," ");  int n = atoi(word);
      if(n > 0) {
        readData(F,n+1,cdf);
        for(int jg=0; jg <= n; jg++) {
          adum[iad] = cdf[jg];
          iad++;
          assert(iad < ltot);
        }
      }
    }
  }
  if(lol[2] > 0) {      // inelastic anisotropy
    if(prlev) printf(" read iad=%d at aniso\n",iad);
    fgets(buf,90,F);
    if(prlev) printf(" at read inelastic aniso : %s",buf);
    int nsg = la[0];
    double pdf[100];
    for(int ig=la[2]-1; ig < la[0]; ig++) {
      fgets(buf,90,F);
      word = strtok(buf," ");  int igg = atoi(word);
      word = strtok(NULL," "); int js = atoi(word);
      assert(igg-1 == ig);
      assert(js > 0);
      readData(F,js,pdf);
      for(int j=0; j < js; j++) {
        adum[iad] = pdf[j];
        iad++;
      }
    }
    if(prlev) printf(" read iad=%d at end aniso\n",iad);
  }
  if(lol[1] > 0) {      // (n,2n) down scatter cdf  and photon n2n cdf
    if(prlev) printf(" read iad=%d at n2n cdf\n",iad);
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    int ndat = ld[1];
    double cdf[100];
    for(int ig=0; ig < la[1]; ig++) {
      fgets(buf,90,F);
      word = strtok(buf," ");   int igr = atoi(word);
      word = strtok(NULL," ");  int ndat2 = atoi(word);
      assert(ndat2 == ld[1]);
      if(ndat2 > 0) {
        readData(F,ndat2+1,cdf);
        for(int jg=0; jg <= ndat2; jg++) {
          adum[iad] = cdf[jg];
          iad++;
        }
      }
    }
    if(prlev) printf(" read iad=%d at end of n2n cdf\n",iad);
  }
  if(lol[4] > 0) {	// (n,3n) downscatter cdf
    if(prlev) printf(" read iad=%d at n3n cdf\n",iad);
    fgets(buf,90,F);
    if(prlev) printf("  %s",buf);
    if(ld[4] > 0) { // to handle exception of 74-W data
      double cdf[100];
      for(int ig=0; ig < la[4]; ig++) {
        fgets(buf,90,F);
        word = strtok(buf," ");   int igr = atoi(word);
        word = strtok(NULL," ");  int ndat = atoi(word);
        assert(ndat == ld[4]);
        if(ndat > 0) {
          readData(F,ndat+1,cdf);
          for(int jg=0; jg <= ndat; jg++) {
            adum[iad] = cdf[jg];
            iad++;
          }
        }
      }
    }
    else { // special handling 
      for(int ig=0; ig < la[4]; ig++) {
        adum[iad] = 1.0;    iad++;
      }
    }
    if(prlev) printf(" read iad=%d at end of n3n cdf\n",iad);

  }
  /*  elastic eq area xmu  */
  if(prlev) printf(" read iad=%d at xmu\n",iad);
  fgets(buf,90,F);
  if(prlev) printf("  %s",buf);
  if(prlev) printf(" iad=%d at eq.xmu\n",iad);
  double xmu[100];
  for(int ig=0; ig < ngrpx; ig++) {
    fgets(buf,90,F);
    word = strtok(buf," ");   int igr = atoi(word);
    word = strtok(NULL," ");  int n = atoi(word);
    adum[iad] = n;	iad++;
    if(n > 0) {
      readData(F,n,xmu);
      for(int jg=0; jg < n; jg++) {
        adum[iad] = xmu[jg]; iad++;
      }
    }
  }
  return iad;
};	// XSiso::microblk


void XSiso::asm_fast(FILE *F, int nfast_gps)
{
  char buf[90];
  char *word;
  int nidar[8];
  fgets(buf,90,F);
  for(int j=0; j < 8; j++) {
    if(j == 0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    nidar[j] = atoi(word);
  }
  if(prlev > 1) {
  printf("nidar:");
  for(int j=0; j < 8; j++) printf(" %d",nidar[j]);
  printf("\n");
  }
  fgets(buf,90,F);
  word = strtok(buf," ");   int ltot = atoi(word);
  word = strtok(NULL," ");  int iwa  = atoi(word);
  word = strtok(NULL," ");  int iwf = atoi(word);
  word = strtok(NULL," ");  int iws = atoi(word);
  word = strtok(NULL," ");  int jframe = atoi(word);
  word = strtok(NULL," ");  int ltype = atoi(word);
  word = strtok(NULL," ");  int iwr = atoi(word);
  if(prlev) printf("ltot=%d iwa=%d iwf=%d iws=%d lframe=%d ltype=%d iwr=%d\n",
        ltot,iwa,iwf,iws,jframe,ltype,iwr);
  if(iwf > 0) {
    printf("*** NO fission allowed for Sera calc. iwf=%d\n",iwf);
    exit(0);
  }
  int lol[5];
  fgets(buf,90,F);
  for(int j=0; j < 5; j++) {
    if(j==0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    lol[j] = atoi(word);
  }
  int la[5];
  fgets(buf,90,F);
  for(int j=0; j < 5; j++) {
    if(j==0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    la[j] = atoi(word);
  }
  int ld[5];
  fgets(buf,90,F);
  for(int j=0; j < 5; j++) {
    if(j==0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    ld[j] = atoi(word);
  }
  if(prlev) {
  printf("lol:");
  for(int j=0; j < 5; j++) printf(" %d",lol[j]);
  printf("\n");
  printf("la :");
  for(int j=0; j < 5; j++) printf(" %d",la [j]);
  printf("\n");
  printf("ld :");
  for(int j=0; j < 5; j++) printf(" %d",ld [j]);
  printf("\n");
  }
  fgets(buf,90,F);
  word = strtok(buf," ");   nrr = atoi(word);
  word = strtok(NULL," ");  awr = atof(word);
  word = strtok(NULL," ");  spot = atof(word);
  if(prlev) printf(" nrr=%d awr=%.3lf spot=%.2le\n",nrr,awr,spot);
  adum_fast = new double[ltot];
  int ltotal = microblk(F,nfast_gps,adum_fast,ltot,iwa,iwf,iws,lol,la,ld);
  assert(ltotal == ltot);
  /*  save for later use : collect  */
  ltot_fast = ltot;
  iwa_fast = iwa;
  iwf_fast = iwf;
  iws_fast = iws;
  for(int i=0; i < 5; i++) lol_fast[i] = lol[i];
  for(int i=0; i < 5; i++) la_fast[i] = la[i];
  for(int i=0; i < 5; i++) ld_fast[i] = ld[i];
};	//  XSiso::asm_fast

void XSiso::asm_thermal(FILE *F, int ntherm_gps)
{
  char buf[90];
  fgets(buf,90,F);
  if(prlev) printf("  therm: %s",buf);
  siga = new double[ntherm_gps];
  readData(F, ntherm_gps,siga);
  fgets(buf,90,F);
  if(prlev) printf("  therm: %s",buf);
  sigs = new double[ntherm_gps];
  readData(F, ntherm_gps,sigs);
  fgets(buf,90,F);
  if(prlev) printf("  therm: %s",buf);
  P0 = new double[ntherm_gps*ntherm_gps];
  for(int jg=0; jg < ntherm_gps; jg++) {
    readData(F,ntherm_gps, P0+jg*ntherm_gps);
  }
  fgets(buf,90,F);
  if(prlev) printf("  therm: %s",buf);
  P1 = new double[ntherm_gps*ntherm_gps];
  for(int jg=0; jg < ntherm_gps; jg++) {
    readData(F,ntherm_gps, P1+jg*ntherm_gps);
  }
};

void XSiso::asm_gamma(FILE *F, int ngam_gps)
{
  char buf[90];
  char *word;
  fgets(buf,90,F);
  if(prlev) printf(" Gamma  %s",buf);
  int nidar[8];
  fgets(buf,90,F);
  for(int i=0; i < 8; i++) {
    if(i==0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    nidar[i] = atof(word);
  }
  if(prlev > 1) {
  printf("G nidar:");
  for(int i=0; i < 8; i++) printf(" %d",nidar[i]);
  printf("\n");
  }
  fgets(buf,90,F);
  word = strtok(buf," ");	int ltot = atoi(word);
  word = strtok(NULL," ");	int iwa = atoi(word);
  word = strtok(NULL," ");      int iwf = atoi(word);
  word = strtok(NULL," ");      int iws = atoi(word);
  word = strtok(NULL," ");      int jframe = atoi(word);
  word = strtok(NULL," ");      int ltype = atoi(word);
  word = strtok(NULL," ");      int iwr = atoi(word);
  if(prlev) printf("G ltot=%d lwa=%d iwf=%d iws=%d jframe=%d ltype=%d iwr=%d\n",
    ltot,iwa,iwf,iws,jframe,ltype,iwr);
  int lol[8];
  fgets(buf,90,F);
  for(int j=0; j < 5; j++) {
    if(j==0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    lol[j] = atoi(word);
  }
  int la[5];
  fgets(buf,90,F);
  for(int j=0; j < 5; j++) {
    if(j==0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    la[j] = atoi(word);
  }
  int ld[5];
  fgets(buf,90,F);
  for(int j=0; j < 5; j++) {
    if(j==0) word = strtok(buf," ");
    else word = strtok(NULL," ");
    ld[j] = atoi(word);
  }
  if(prlev) {
  printf("G lol:");
  for(int j=0; j < 5; j++) printf(" %d",lol[j]);
  printf("\n");
  printf("G la :");
  for(int j=0; j < 5; j++) printf(" %d",la [j]);
  printf("\n");
  printf("G ld :");
  for(int j=0; j < 5; j++) printf(" %d",ld [j]);
  printf("\n");
  }
  /*  we dont take awr from gamma  */
  fgets(buf,90,F);
  word = strtok(buf," ");   nrr = atoi(word);
  word = strtok(NULL," ");  double awr = atof(word);
  word = strtok(NULL," ");  double spot = atof(word);
  if(prlev) printf(" nrr=%d awr=%.3lf spot=%.2le\n",nrr,awr,spot);
  adum_gam = new double[ltot];
  int ltotal = microblk(F,ngam_gps,adum_gam,ltot,iwa,iwf,iws,lol,la,ld);
  assert(ltotal == ltot);
  /*  save for later use : collect  */
  ltot_gam = ltot;
  iwa_gam = iwa;
  iwf_gam = iwf;
  iws_gam = iws;
  for(int i=0; i < 5; i++) lol_gam[i] = lol[i];
  for(int i=0; i < 5; i++) la_gam[i] = la[i];
  for(int i=0; i < 5; i++) ld_gam[i] = ld[i];


};	//  XSiso::asm_gamma
