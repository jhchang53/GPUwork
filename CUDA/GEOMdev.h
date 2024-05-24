/*
	GEOMdev.h
*/
#ifndef _inc_Ray
#define _inc_Ray
struct Ray {
  int here,next;
  double dist;
};
#endif

__device__ void TRAN_set(double target[], double Zb, double theta, double phi, double psi);

__device__ Ray howfar(double P[]);
