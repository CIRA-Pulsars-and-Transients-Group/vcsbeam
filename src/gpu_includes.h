#ifndef _GPU_INCLUDES_H__
#define _GPU_INCLUDES_H__

// #ifndef __NVCC__
// #define __NVCC__ // should be set by compiler !!!
// #endif


#ifdef __NVCC__
   #include <cuda_runtime.h>
   #include <cufft.h>
   #include <cuComplex.h>   
#else
   #include <hipfft/hipfft.h>
#endif

#endif
