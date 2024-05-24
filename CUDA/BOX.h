/*
	BOX.h
*/
__device__ void BOX_set(double bmin[], double bmax[]);
__device__ int BOX_DoRayQuery(const double Pin[], double t[2]);

