#if defined(__XSLdev)
__device__ int num_fast_gps,num_therm_gps,num_neut_gps,num_gam_gps,num_iso;
__device__ double *ev_neut,*ev_gam;

__device__ int *id_fast,*id_therm,*id_gam, *jscat,*nbound;
__device__ double *awr;		// atomic mass
__device__ int *ind_fast;	// index to fast_cdf
__device__ double *fast_cdf;	// band 5 reaction cdf array

__device__ int *ind_xmu;	// index to xmu size: niso*fast_gps;
__device__ double *xmu;

__device__ int *ind_inel,*inel_grps,*inel_cdfs;
__device__ double *cdf_inel,*aniso_inel;

__device__ int *ind_n2n,*n2n_grps,*n2n_cdfs;
__device__ double *cdf_n2n;

__device__ int *ind_n3n,*n3n_grps,*n3n_cdfs;
__device__ double *cdf_n3n;

__device__ int *ind_therm;
__device__ double *therm_cdf;	//  band 2 reaction cdf array (capt+elas)

__device__ int *ind_gam;
__device__ double *gam_cdf;	// band 3 reaction type cdf

__device__ double *sig_2200,*sig_gprod;
__device__ int *ind_gprod;
__device__ int *num_gamline;
__device__ double *g_energy,*g_yield;

__device__ double *xkerma_n,*xkerma_g;

__device__ void checkFast(int iso);
__device__ void checkXmu(int iso);
__device__ void checkInel(int iso);
__device__ void checkN2n(int iso);
__device__ void checkN3n(int iso);
__device__ void checkTherm(int iso);
__device__ void checkGam(int iso);
__device__ void checkGprod(int iso);
__device__ void checkKerma(int iso);
#else
extern __device__ int num_fast_gps,num_therm_gps,num_neut_gps,num_gam_gps,num_iso;
extern __device__ double *ev_neut,*ev_gam;
extern __device__ int *id_fast,*id_therm,*id_gam, *jscat,*nbound;
extern __device__ double *awr;
extern __device__ int *ind_fast;       // index to fast_cdf
extern __device__ double *fast_cdf;    // band 5 reaction cdf array

extern __device__ int *ind_xmu;        // index to xmu
extern __device__ double *xmu;

extern __device__ int *ind_inel,*inel_grps,*inel_cdfs;
extern __device__ double *cdf_inel,*aniso_inel;

extern __device__ int *ind_n2n,*n2n_grps,*n2n_cdfs;
extern __device__ double *cdf_n2n;

extern __device__ int *ind_n3n,*n3n_grps,*n3n_cdfs;
extern __device__ double *cdf_n3n;

extern __device__ int *ind_therm;
extern __device__ double *therm_cdf;

extern __device__ int *ind_gam;
extern __device__ double *gam_cdf;   // band 3 reaction type cdf

extern __device__ double *sig_2200,*sig_gprod;
extern __device__ int *ind_gprod;
extern __device__ int *num_gamline;
extern __device__ double *g_energy,*g_yield;

extern __device__ double *xkerma_n,*xkerma_g;


__global__ void XSLset(int look);
#endif
