#ifndef __GPUFFT_H__
#define __GPUFFT_H__

#ifdef __NVCC__
  // NVIDIA / CUFFT :

  #include <cufftw.h>
  #include <cufft.h>
  #define gpufftResult   cufftResult
  #define gpufftDestroy  cufftDestroy
  #define GPUFFT_SUCCESS CUFFT_SUCCESS
  #define gpufftComplex  cufftComplex
  #define gpufftPlanMany cufftPlanMany
  #define gpufftHandle   cufftHandle
  #define gpufftPlan2d   cufftPlan2d
  #define gpufftExecC2C  cufftExecC2C
  #define GPUFFT_C2C     CUFFT_C2C
  #define GPUFFT_FORWARD CUFFT_FORWARD
#else
  // AMD / HIP :
  
  #include <hipfft/hipfft.h>
  
  #define gpufftResult   hipfftResult
  #define gpufftDestroy  hipfftDestroy
  #define GPUFFT_SUCCESS HIPFFT_SUCCESS
  #define gpufftComplex  hipfftComplex
  #define gpufftPlanMany hipfftPlanMany
  #define gpufftHandle   hipfftHandle
  #define gpufftPlan2d   hipfftPlan2d
  #define gpufftExecC2C  hipfftExecC2C
  #define GPUFFT_C2C     HIPFFT_C2C
  #define GPUFFT_FORWARD HIPFFT_FORWARD
#endif


#endif
