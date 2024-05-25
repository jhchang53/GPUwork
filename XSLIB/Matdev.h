#if defined(__Matdev)
__device__ int nmater;
__device__ int *index_iso,*index_fmac,*index_fcdf,*index_tmac,*index_tcdf;
__device__ int *iso;
__device__ double *fast_sigtot,*fast_cdf,*therm_sigtot,*therm_cdf;

#else
extern __device__ int nmater;
extern __device__ int *index_iso,*index_fmac,*index_fcdf,*index_tmac,*index_tcdf;
extern __device__ int *iso;
extern __device__ double *fast_sigtot,*fast_cdf,*therm_sigtot,*therm_cdf;

__global__ void MatSet();
#endif
