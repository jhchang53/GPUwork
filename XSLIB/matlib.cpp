/*
        matlib.cpp
        read material file
	and create mater sigtot and cdf
*/
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "XSLIB.h"
#include "Mater.h"

using namespace std;

int main()
{
  XSLIB *xslib = new XSLIB();
  XSL *xsl = xslib->openlib("seraMC.libupd4x.punch");

  Mater *mat = new Mater();
  mat->readMat("added2.mat");
  mat->setXSL(xsl);
  vector<string> matlist = {"buffer","tissue","PMMA","brain_ICRU_46_adult"};
  mat->setupMatLib(matlist);
  // mat->print(2);
  mat->sendMat();
  delete mat;
  delete xslib;
};
