/*
	mater.cpp
	read material file
*/
#include <stdio.h>
#include <stdlib.h>
#include "XSLIB.h"
#include "Mater.h"

int main()
{
  XSLIB *xslib = new XSLIB();
  XSL *xsl = xslib->openlib("seraMC.libupd4x.punch");

  Mater *mat = new Mater();
  mat->readMat("added2.mat");
  mat->setXSL(xsl);
  int m = mat->select("brain_ICRU_46_adult");
  // int m = mat->select("tissue");
  // int m = mat->select("buffer");
  // int m = mat->select("void_true");
  mat->collect(m);
  delete mat;
  delete xslib;
};
