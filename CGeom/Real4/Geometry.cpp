/*
	Geometry.cpp - float version
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Geometry.h"

Geometry::Geometry()
{
  prlev = 0;
  trc = nullptr;
  inUV = 0;
  nextzone = 1;
  ncgs = 0;
  INF = 1.0/0.0;
  isok = nullptr;
  t_min = nullptr;
  t_max = nullptr;
  chx = ' ';
  box = nullptr;
  vox = nullptr;
};

Geometry::~Geometry()
{
  delete [] isok;
  delete [] t_min;
  delete [] t_max;

  delete box;
  delete vox;
};

void Geometry::setPrlev(int prl)
{
  prlev = prl;
};	// Geomtry::setPrlev

void Geometry::setTR(TRcoord *trc_i)
{
  trc = trc_i;
};

Ray Geometry::track(const float p[], const float w[])
{
  Ray rc;
  /*  p, w is UV  */
  char chx = vox->getChar(p);
  inUV  = 0;
  if(chx >= 32) inUV = 1;
  if(inUV) {
    if(prlev) printf("== %s:%d inside of UV\n",__FILE__,__LINE__);
    if(chx < 32) chx = 32;	
    rc = vox->track(p,w,chx);
    if(prlev) {
      printf("  Vox here=%c next=%c dist=%.2le\n",chx,rc.next, rc.dist);
    }
    chx = rc.next;
    if(chx < 0) {
      inUV = 0;
      return rc;
    }
    // printf(" at line %d  chx=%c dist=%.2le\n",__LINE__,chx,rc.dist);
    if(chx != ' ')  return rc;
    // printf("  passing line %d\n",__LINE__);
    /*  next cell is void  */
    if(prlev) printf("== %s:%d next cell is void\n",__FILE__,__LINE__);
    /*  check all CG if there is intersection */
    trc->U2S(p,CGp);
    trc->U2Svec(w,CGw);
    Ray cgray = CGtrack(CGp,CGw);
    if(prlev) printf("CG  ray=%d %d %.3lf\n",cgray.here,cgray.next,cgray.dist);
    // int isout = vox->checkXYZ();
    if(cgray.here == 0) return rc;
    if(prlev) printf("== %s:%d next cell is void\n",__FILE__,__LINE__);
    rc.next = -1;
    if(prlev) printf("**  escape\n");
    inUV = 0;
    // printf("  passing line %d\n",__LINE__);
    rc.dist = 0.0;
    return rc;
  }
  else {	// not in a UV box
    /* it can be CG or UV  */
    // printf("  passing line %d\n",__LINE__);
    trc->U2S(p,CGp);
    trc->U2Svec(w,CGw);
    Ray cgray = CGtrack(CGp,CGw);
    if((cgray.here == 0) && (cgray.dist > 1.0e+10)) {
      // printf("  passing line %d\n",__LINE__);
      // check UV box */
      float t[2];
      int meetBox = box->DoRayQuery(p,w, t);
      if(prlev) printf("Box  meet=%d t0=%.3lf t1=%.3lf\n",meetBox,t[0],t[1]);
      if(!meetBox) {
        rc.here = -1;  rc.next = -1;
        rc.dist = 0.0;
        if(prlev) printf("**  escape\n");
        return rc;
      }
      else {
        rc.here = 0;
        rc.next = 0;
        rc.dist = t[0];
        if(rc.dist < 0.0) rc.dist = 0;
        inUV = 1;
        if(prlev) printf("#  entering Box at %.2le\n",t[0]);
        return rc;
      }
    }
    else {
      if(cgray.dist < 1.0e-10) cgray.dist = 1.0e-10;
      return cgray;
    }
  }
  printf("== %s:%d\n",__FILE__,__LINE__);
  exit(0);
};

void Geometry::setCG(int ncgs_i, CG **cg_i)
{
  ncgs = ncgs_i;
  cg = cg_i;
  /* reserve spaces */
  isok = new int[ncgs];
  t_min = new float[ncgs];
  t_max = new float[ncgs];
};

Ray Geometry::CGtrack(float p[], float w[])
{
  Ray ray;
  /*  find track lengths of all CG */
  int thiszone = 0;
  float dist = 1.0e+20;
  int zprod = 1;
  float trk[2];
  for(int n=0; n < ncgs; n++) {
    CG *cgeom = cg[n];
    int ok = cgeom->track(p,w, trk);
    if(prlev >1) printf("CG%d: ok=%d trk=%.3e %.3e\n",n,ok,trk[0],trk[1]);
    isok[n] = ok;
    if(ok) {
      t_min[n] = trk[0];
      t_max[n] = trk[1];
      /*  find this zone */
      if((trk[0] <= 0.0) && (trk[1] > 0.0)) {
        thiszone += zprod;
      }
      if((t_min[n] > 0.0) && (dist > t_min[n])) dist = t_min[n];
      if((t_max[n] > 0.0) && (dist > t_max[n])) dist = t_max[n];
    }
    zprod = 2*zprod;	// prepare for next CG
  }
  /*  find next zone  */
  if(prlev) printf(" CGtrack  thiszone=%d dist=%.3le\n",thiszone,dist);
  /*  check next zone  */
  zprod = 1;
  int nextzone = 0;
  for(int n=0; n < ncgs; n++) {
    if(isok[n]) {
      if((t_min[n] <= dist) && (dist < t_max[n])) {
        nextzone += zprod;
      }
    }
    zprod = 2*zprod;
  }
  if(prlev) printf(" CGtrack  nextzone=%d\n",nextzone);
  ray.here = thiszone;
  ray.next = nextzone;
  ray.dist = dist;
  return ray;
};	// Geometry::CGtrack

/*  UV geometry  */
void Geometry::setLattice(int nxyz[], float dxyz[], char *uv)
{
  float bmin[] = {0.0,0.0,0.0};
  float bmax[3];
  for(int i=0; i < 3; i++) {
    bmax[i] = nxyz[i]*dxyz[i];
  }
  box = new Box();
  box->set(bmin,bmax);
  vox = new Vox();
  vox->setLattice(nxyz,dxyz,uv);
};	// Geometry::setLattice

Ray Geometry::UVtrack(float p[], float w[], char chx)
{
  Ray ray = vox->track(p,w,chx);
  return ray;
};	// Geometry::UVtrack
