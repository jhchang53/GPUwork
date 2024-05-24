/*
	Scatter.cpp
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Scatter.h"
#include "Xor128.h"

Scatter::Scatter()
{
  prlev = 0;
  PI = 2*acos(0.0);
};

Scatter::~Scatter()
{
};

void Scatter::setPrlev(int prl)
{
  prlev = prl;
};	// Scatter::setPrlev

double Scatter::isotropic(const double A, const double E, double w[])
{
  /* isotropic in CMS  */
  double cms = 1.0-2.0*ranf();
  /* compute cosine at lab system  */
  double xmu = (A*cms+1.0)/sqrt(A*(A+2*cms)+1.0);
  if(prlev) {
    printf("cms=%.3lf ",cms);
    printf(" xmu=%.4lf\n",xmu);
  }
  /*  rotate  */
  double wvec[3];
  rotate(xmu,w,wvec);
  for(int i=0; i < 3; i++) w[i] = wvec[i];
  /*  compute energy after scattering  */
  double ratio = (xmu+sqrt(xmu*xmu+A*A-1.0))/(A+1.0);
  return E*ratio*ratio;
};	// Scatter::isotropic

double Scatter::P1(const double p0, const double p1, const double w[], double wvec[])
{
  /*  P0, P1 distribution  */
  /*  use rejection method to avoid sqrt */
  double a = 1-p1/p0;
  double b = p1/p0;
  double A = a+fabs(b);
  double u,x;
  do {
    u = ranf();
    x = ranf();
  }
  while (u*A > a+b*x);
  double xmu = 2*x-1;
  if(prlev) printf(" xmu=%.5lf\n",xmu);
  return xmu;
};	// Scatter::P1

double Scatter::Compton(const double ein, double &enew)
{
  double x,xmu;
  double enew;
  double ra1,ra2,ra3;
  double alpha = ein/0.511034e+6;

  double t4 = 2*alpha;
  double t1 = 1.0/alpha;
  double t2,t3;
  if(ein <= 1.5e+6) {
    /*  Kahn's method  */
    /* ref. Herman Kahn, Application of Monte Carlo, RM-1237-AEC,  p.64 (1965)  */
    t3 = (t4+1.0)/(t4+9.0);

    while(1) {
      ra1 = ranf();
      ra2 = ranf();
      ra3 = ranf();
      if(ra1 > t3) {
        x = 1.0/(1.0+t4*ra2);
        t2 = 4.0*(x-x*x);
        printf(" %s:%d t2=%.2le\n",__FILE__,__LINE__,t2);
        if(ra3 <= t2) {
          xmu = t1 + 1.0-t1/x;
          printf(" %s:%d xmu=%.2le\n",__FILE__,__LINE__,xmu);
          break;
        }
      }
      else {
        x = (1.0+t4*ra2)/(t4+1.0);
        xmu = t1 + 1.0 - t1/x;
        t2 = 0.5*(xmu*xmu+x);
        printf(" %s:%d t2=%.2le\n",__FILE__,__LINE__,t2);
        if(ra3 <= t2) break;
      }
    }
    enew = ein*x;
    return xmu;
  }
  else {
    /*  Koblinger's method  */
#ifdef XX
    ra1 = ranf();
    ra2 = ranf();
    double a = t1*t1;
    double beta = 1.0 + t4;
    t3 = 1.0/(beta*beta);
    double h = 1.0/(4.0*t1+alpha*(1.0+beta)*t3 - a*(beta+1.0-alpha*alpha)*log(beta));
    t2 = (t1+t1)/h;
    if(ra1 <= t2) {
      x = 1.0/(1.0+t4*ra2);
    }
    else {
      t2 = 2*t2;
      if(ra1 <= t2) {
        x = (1.0+t4*ra2)/beta;
      }
      else {
        double gamma = 1.0-t3;
        t2 = t2 + 0.5*h*gamma;
        if(ra1 <= t2) {
          x = sqrt(1.0-gamma*ra2);
        }
        else {
          x = pow(beta,-ra2);
        }
      }
    }
xmu = t1 + 1.0 - t1/x;
    xmu = t1 + 1.0 -t1/4;
#else
    /*  L.Koblinger, Direct Samplin from the Klein-Nishina Didstribution */
    /*  for Photon Energies Above 1.4 MeV, NSE 56(2) pp.218-219, 1979  */
    /*  use prob. mixing method, Lux p.9  */
    double alpha2 = alpha*alpha;
    double A = 1/alpha2;
    double B = 1-2*(alpha+1)/alpha2;
    double C = (1+2*alpha)/alphasq;
    double D = 1;
    double p1 = 2/alpha;
    double beta = 1+2*alpha;
    double p2 = (1-(1+beta))*log(beta);
    double p3 = p1;
    double gamma = 1-1.0/(beta*beta);
    double p4 = gamma/2;
    double H = 1.0/(p1+p2+p3+p4);

    double s = ranf();
    double x;
    if(s < H*p1) {
      x = 1+2*alpha*ranf();
    }
    else if(s < H*(p1+p2)) {
      x = pow(beta,ranf());
    }
    else if(s < H*(p1+p2+p3)) {
      x = beta/(1+2*alpha*ranf());
    }
    else {
      x = 1.0/sqrt(1-gamma*ranf());
    }
    double eout = ein/x;
    double xmu = t1 + 1.0 - t1/x;
    return xmu;
    
#endif
  }
};	//  Scatter::Compton

void Scatter::rotate(const double xmu, double k[], double vrot[])
{
  /*  find rotation axe  */
  int jaxe=0;
  double amp = 1.0;
  for(int i=0; i < 3; i++) {
    if(fabs(k[i]) < amp) {
     jaxe = i;
     amp = fabs(k[i]);
    }
  }
  if(prlev) {
    printf(" w=(%.4lf,%.4lf,%.4lf)\n",k[0],k[1],k[2]);
    printf(" jaxe=%d\n",jaxe);
  }
  /*  find a vector normal to k and jaxe then rotate around jaxe angle xmu  */
  double c = xmu;
  double s = sqrt(1.0-c*c);
  double v[3];
  if(jaxe == 0) {
    double knorm = sqrt(k[1]*k[1]+k[2]*k[2]);
    v[0] =  k[0]*c - knorm*s;
    v[1] =  k[1]*c + k[0]*k[1]*s/knorm;
    v[2] =  k[2]*c + k[2]*k[0]*s/knorm;
  }
  else if(jaxe == 1) {
    double knorm = sqrt(k[2]*k[2]+k[0]*k[0]);
    v[1] =  k[1]*c - knorm*s;
    v[2] =  k[2]*c + k[1]*k[2]*s/knorm;
    v[0] =  k[0]*c + k[0]*k[1]*s/knorm;
  }
  else if(jaxe == 2) {
    double knorm = sqrt(k[0]*k[0]+k[1]*k[1]);
    v[2] =  k[2]*c - knorm*s;
    v[0] =  k[0]*c + k[2]*k[0]*s/knorm;
    v[1] =  k[1]*c + k[1]*k[2]*s/knorm;
  }
  else abort();
  if(prlev) {
    double kk = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    double kv = k[0]*v[0] + k[1]*v[1] + k[2]*v[2];
    double vv = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    printf("  c=%.5lf s=%.5lf",c,s);
    printf("  kk=%.5lf vv=%.5lf kv=%.5lf\n",kk,vv,kv);
    printf(" v=(%.4lf,%.4lf,%.4lf)\n",v[0],v[1],v[2]);
  }
  /*  use Rodrigue's formula */
  double dotkv = k[0]*v[0] + k[1]*v[1] + k[2]*v[2];
  double kvprod[3];
  kvprod[0] = k[1]*v[2]-k[2]*v[1];
  kvprod[1] = k[2]*v[0]-k[0]*v[2];
  kvprod[2] = k[0]*v[1]-k[1]*v[0];

#ifdef SIMPLE
  double theta = 2*PI*ranf();
  double costh = cos(theta);
  double sinth = sin(theta);
#else
  /*  Neumann method  Lux p.21*/
  double rho1,rho2,rsq;
  while(1) {
    rho1 = 2*ranf()-1.0;
    rho2 = 2*rand()-1.0;
    rsq = rho1*rho1 + rho2*rho2;
    if(rsq <= 1.0) break;
  }
  double rnorm = sqrt(rsq);
  double costh = rho1/rnorm;
  double sinth = rho2/rnorm;
#endif
  for(int i=0; i < 3; i++) {
    vrot[i] = v[i]*costh + kvprod[i]*sinth + k[i]*dotkv*(1.0-costh);
  }
  if(prlev) {
    /*  check result  */
    printf(" vrot=(%.4lf,%.4lf,%.4lf)\n",vrot[0],vrot[1],vrot[2]);
    double klen = k[0]*k[0]+k[1]*k[1]+k[2]*k[2];
    double vlen = vrot[0]*vrot[0]+vrot[1]*vrot[1]+vrot[2]*vrot[2];
    printf("  klen=%.5lf", klen);
    printf("  vlen=%.5lf", vlen);
    double kvrot = k[0]*vrot[0] + k[1]*vrot[1] + k[2]*vrot[2];
    printf("  kvrot=%.5lf",kvrot);
    printf("  err=%.3le\n",kvrot-xmu);
  }
};	// Scatter::rotate

