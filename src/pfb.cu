/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>

#include <mwalib.h>

#include "vcsbeam.h"

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    /* Wrapper function for GPU/CUDA error handling. Every CUDA call goes through
       this function. It will return a message giving your the error string,
       file name and line of the error. Aborts on error. */

    if (code != 0)
    {
        fprintf(stderr, "GPUAssert:: %s - %s (%d)\n", cudaGetErrorString(code), file, line);
        if (abort)
        {
            exit(code);
        }
    }
}

// define a macro for accessing gpuAssert
#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__, true);}

/*******************************
 * FORWARD (ANALYSIS) FINE PFB *
 *******************************/

/* The following kernels are part of an "offline" version of the "fine PFB"
   algorithm that was implemented on the FPGAs of Phase 1 & 2 of the MWA. As
   described in McSweeney et al. (2020), this algorithm is a version of the
   "weighted overlap-add" algorithm (see their Eq. (3)):

              K-1
     X_k[m] = SUM b_m[n] e^(-2πjkn/K),
              n=0

   where

              P-1
     b_m[n] = SUM h[Kρ-n] x[n+mM-Kρ].
              ρ=0

   A full description of these symbols is given in the reference, but it
   should be noted here that, for the kernels below,

     - "x" represents the input data (INDATA),
     - "X" represents the output data (OUTDATA),
     - "h" represents the filter coefficients (FILTER_COEFFS),
     - "P" represents the number of taps (NTAPS)
     - "b" represents the weighted overlap-add array (WEIGHTED_OVERLAP_ARRAY)

   The algorithm is broken up into three parts:
     (1) the weighted overlap-add, which forms "b" from "x" and "h",
     (2) the FFT, which is implemented using the cuFFT library, and
     (3) the demotion and packaging of the result into the required
         output format, which is described in the Appendix of McSweeney et al.
         (2020).

   Notes:

    - INDATA is expected to have the same data layout as delivered by mwalib's
      mwalib_voltage_context_read_second() function when applied to MWAX high
      time resolution data, whose format description can be found at
      https://wiki.mwatelescope.org/display/MP/MWA+High+Time+Resolution+Voltage+Capture+System
      The samples are organised according to the vMWAX_IDX macro.

    - OUTDATA will have the same data layout as "recombined" legacy VCS data:
      (4+4)-bit complex samples, with imaginary component occupying the first
      4 bits, and samples organised according to the v_IDX macro.

    - Each thread will operate on one RF input (i.e. antenna/pol combination)
      and generate the spectrum for a single "fine-channelised" time step
      (the "m" index). The "weighted overlap-add" array ("b") also resides in
      device memory, in order that it can make use of cuFFT. This is a 
 */

__global__ void legacy_pfb_weighted_overlap_add( char2 *indata,
        int *filter_coeffs, void *weighted_overlap_array )
{
    // Parse the block and thread idxs:
    //   <<<(nspectra,I,P),(K)>>>
    // where nspectra is the number of output spectra,
    //       I        is the number of RF inputs,
    //       K        is the size of the output spectrum
    //       P        is the number of taps
    // ...and put everything into the mathematical notation used in McSweeney
    // et al. (2020) (see equation in comments above)
    //int nspectra = gridDim.x;
    int    I = gridDim.y;
    int    P = gridDim.z;

    int    K = blockDim.x;

    int    m = blockIdx.x;
    int    i = blockIdx.y;
    int    p = blockIdx.z;

    int    M = K; // This enforces a critical sampled PFB
    int    n = threadIdx.x;

    int   *h = filter_coeffs;
    char2 *x = indata;
    int2  *b = (int2 *)weighted_overlap_array;

    // Because we have included a whole voltage block at the beginning of the
    // buffer, index "m" actually needs to start at N-P, where N is the number
    // of M-length strides in a voltage block, and P is the number of taps.
    // For MWAX data, N is always 500.
    int mprime = m + 500 - P;

    // Now calculate the index into the various arrays that also
    // takes into account the fact that these arrays contain all RF inputs.
    // MEMO TO SELF: My current going theory is that I don't have to do any
    // re-ordering of the antennas, as that is dealt with elsewhere.
    int h_idx = K*P - K*p - n - 1; // This reverses the filter (think "convolution")
    unsigned int x_idx = vMWAX_IDX((unsigned int)(mprime*M + p*K + n), i, I);
    unsigned int b_idx = m*(K*I) + i*(K) + n; // This puts each set of K samples to
                                              // be FFT'd in a contiguous memory block

    // Now perform the weighted overlap-add operation
    int   hval = h[h_idx];
    char2 xval = x[x_idx];

    // In keeping with the original FPGA implementation, the result now needs to
    // be demoted and rounded.
    int X = hval*(int)xval.x;
    int Y = hval*(int)xval.y;

    // Promote the result to floats and put it in the b array in global memory
    // in preparation for being FFTed
    atomicAdd( &(b[b_idx].x), X );
    atomicAdd( &(b[b_idx].y), Y );

//if (m == 0 && i == 0) printf( "%u %d %d %d %d %d %d %d %d %d %d %d %d\n", n, p, K, P, h_idx, hval, xval.x, xval.y, X, Y, b_idx, b[b_idx].x, b[b_idx].y );
}

__global__ void fpga_rounding_and_demotion( void *data )
{
    // Parse the block and thread idxs:
    //   <<<(nspectra,I),(K)>>>
    // where nspectra is the number of output spectra,
    //       I        is the number of RF inputs,
    //       K        is the size of the output spectrum
    // ...and put everything into the mathematical notation used in McSweeney
    // et al. (2020) (see equation in comments above)
    //int nspectra = gridDim.x;
    int      I = gridDim.y;
    int      K = blockDim.x;
    int      m = blockIdx.x;
    int      i = blockIdx.y;
    int      n = threadIdx.x;

    int2    *b = (int2 *)data;

    unsigned int b_idx = m*(K*I) + i*(K) + n;

    int X = b[b_idx].x;
    int Y = b[b_idx].y;

    // Rounding:
    if (X > 0)  X += 0x2000;
    if (Y > 0)  Y += 0x2000;

    // Demotion:
    X >>= 14;
    Y >>= 14;

    // Put the result back into global memory as (32-bit) floats
    float *fx = (float *)&(b[b_idx].x);
    float *fy = (float *)&(b[b_idx].y);

    *fx = (float)X;
    *fy = (float)Y;
}

__global__ void int2float( void *data, double scale )
{
    // Parse the block and thread idxs
    // Each thread handles a single int
    // It is assumed that floats are always the same size as ints (=32 bits)
    int    i           = blockIdx.x*blockDim.x + threadIdx.x;
    int   *data_as_int = (int *)data;
    int   *pint        = &(data_as_int[i]);
    float *pfloat      = (float *)pint;

    // Put the result back into global memory as (32-bit) float, with optional scaling factor applied
    *pfloat = (float)(*pint) * scale;
}

__global__ void pack_into_recombined_format( cuFloatComplex *ffted, void *outdata, int *i_idx, pfb_flags flags )
/* This is the final step in the forward fine PFB algorithm that emulates what
   was implemented on the MWA FPGAs in Phase 1 & 2 (see above for details).
   At this point, the FFTED array contains the Fourier-transformed data that
   already represents the final channelisation. All that remains to be done is
   to pack it into the same format as the VCS recombined data.

   The available flags for PACKING_FLAGS are:
       PFB_COMPLEX_INT4    = pack into 4+4-bit complex signed ints
       PFB_COMPLEX_FLOAT64 = pack into 64+64-bit complex signed floats (= cuDoubleComplex)

   Kernel signature:
     <<<(nspectra,K),I>>>
   where
     nspectra is the number of (fine-channelised) time samples
     K        is the number of channels
     I        is the number of RF inputs
*/
{
    // Parse the kernel signature, using the same mathematical notation
    // described above
    //int nspectra = gridDim.x;
    int K        = gridDim.y;
    int m        = blockIdx.x;
    int k        = blockIdx.y;
    int i        = threadIdx.x;
    int I        = blockDim.x;

    int i_out    = i_idx[i];

    int kprime   = (k + K/2) % K; // Puts the DC bin in the "middle"

    cuFloatComplex  *b = ffted;

    // Calculate the idxs into b and X
    unsigned int b_idx = m*(K*I) + i*(K) + k;
    unsigned int X_idx = v_IDX(m, kprime, i_out, K, I);

    // Pull the values to be manipulated into register memory (because the
    // packing macro below involves a lot of repetition of the arguments)
    // The division by K is to normalise the preceding FFT
    double re = b[b_idx].x / K;
    double im = b[b_idx].y / K;

    // Put the packed value back into global memory at the appropriate place
    if (flags & PFB_COMPLEX_INT4)
    {
        uint8_t *X = (uint8_t *)outdata;
        if (flags & PFB_IMAG_PART_FIRST)
            X[X_idx] = PACK_NIBBLES(im, re);
        else
            X[X_idx] = PACK_NIBBLES(re, im);
    }
    else // Currently, default is cuDoubleComplex
    {
        cuDoubleComplex *X = (cuDoubleComplex *)outdata;
        if (flags & PFB_IMAG_PART_FIRST)
            X[X_idx] = make_cuDoubleComplex( im, re );
        else
            X[X_idx] = make_cuDoubleComplex( re, im );
    }

    __syncthreads();
}

void vmInitForwardPFB( vcsbeam_context *vm, int M, int nchunks, pfb_flags flags )
/* Create and initialise a forward_pfb struct.
 * A filter must already be loaded into VM->ANALYSIS_FILTER
 * **WARNING! The filter will be forcibly typecast to int!!**
 *
 * Inputs:
 *
 *   VM                - vcsbeam metadata, derived from mwalib
 *   M                 - the stride of the application of the filter
 *                       (this determines the time resolution of the
 *                       channelised data). Set M = K for critically
 *                       sampled PFB.
 *   NCHUNKS           - Split each second into this many GPU processing chunks
 *   FLAGS             - Various flags for whether to allocate memory
 *                       (see pfb_flags enum)
 *
 */
{
    // Create the struct in memory
    vm->fpfb = (forward_pfb *)malloc( sizeof(forward_pfb) );
    forward_pfb *fpfb = vm->fpfb; // (Shorthand variable)

    // Set the next gps second to read to be the first one
    fpfb->flags           = flags; // It doesn't matter if there is "extra" information
                                   // in these flag bits apart from the output format

    // Some of the data dimensions
    unsigned int nsamples = vm->vcs_metadata->num_samples_per_voltage_block *
                            vm->vcs_metadata->num_voltage_blocks_per_second;

    fpfb->nspectra = nsamples / M;
    fpfb->M        = M;
    fpfb->K        = vm->analysis_filter->nchans;
    fpfb->I        = vm->obs_metadata->num_rf_inputs;
    fpfb->P        = vm->analysis_filter->ncoeffs / fpfb->K;
    // ^^^ The user is responsible for making sure that the number of desired
    // channels divides evenly into the number of filter coefficients. No
    // error or warning is generated otherwise, not even if the inferred
    // number of taps (P) is 0.

    // Set up the idxs for the "rf input" output order,
    // and copy to device
    gpuErrchk(cudaMallocHost( (void **)&(fpfb->i_output_idx),   fpfb->I*sizeof(int) ));
    gpuErrchk(cudaMalloc(     (void **)&(fpfb->d_i_output_idx), fpfb->I*sizeof(int) ));

    int i, mwax_idx, legacy_idx;
    uint32_t ant;
    char pol;
    for (i = 0; i < fpfb->I; i++)
    {
        ant = vm->obs_metadata->rf_inputs[i].ant;
        pol = *(vm->obs_metadata->rf_inputs[i].pol);

        mwax_idx   = 2*ant + (pol - 'X');
        legacy_idx = vm->obs_metadata->rf_inputs[i].vcs_order;

        fpfb->i_output_idx[mwax_idx] = legacy_idx;
    }

    gpuErrchk(cudaMemcpyAsync( fpfb->d_i_output_idx, fpfb->i_output_idx, fpfb->I*sizeof(int), cudaMemcpyHostToDevice ));

    // Work out the sizes of the various arrays
    while (fpfb->nspectra % nchunks != 0)
        nchunks++;
    fpfb->chunks_per_second  = nchunks;
    fpfb->chunk              = 0;
    fpfb->nspectra_per_chunk = fpfb->nspectra / nchunks;

    int num_voltage_blocks_per_chunk = vm->vcs_metadata->num_voltage_blocks_per_second / nchunks;
    fpfb->d_htr_size = (num_voltage_blocks_per_chunk + 1) * vm->vcs_metadata->voltage_block_size_bytes;

    if (flags & PFB_COMPLEX_INT4)
    {
        fpfb->vcs_size = fpfb->nspectra * vm->obs_metadata->num_rf_inputs * fpfb->K * sizeof(uint8_t);
    }
    else // i.e., default is cuDoubleComplex
    {
        fpfb->vcs_size = fpfb->nspectra * vm->obs_metadata->num_rf_inputs * fpfb->K * sizeof(cuDoubleComplex);
    }

    fpfb->d_vcs_size = fpfb->vcs_size / nchunks;

    fpfb->htr_stride = (vm->vcs_metadata->num_voltage_blocks_per_second*vm->vcs_metadata->voltage_block_size_bytes)/nchunks;
    fpfb->vcs_stride = fpfb->d_vcs_size;

    fpfb->char2s_per_second =
        (vm->vcs_metadata->num_voltage_blocks_per_second *
         vm->vcs_metadata->voltage_block_size_bytes) / sizeof(char2);
    fpfb->bytes_per_block = vm->vcs_metadata->voltage_block_size_bytes;

    fpfb->weighted_overlap_add_size = vm->v->buffer_size * (sizeof(cuFloatComplex) / sizeof(char2)) / nchunks;
    size_t filter_size = vm->analysis_filter->ncoeffs * sizeof(int);

    // Allocate memory for filter and copy across the filter coefficients,
    // casting to int
    gpuErrchk(cudaMallocHost( (void **)&(fpfb->filter_coeffs),   filter_size ));
    gpuErrchk(cudaMalloc(     (void **)&(fpfb->d_filter_coeffs), filter_size ));
    for (i = 0; i < vm->analysis_filter->ncoeffs; i++)
        fpfb->filter_coeffs[i] = (int)vm->analysis_filter->coeffs[i]; // **WARNING! Forcible typecast to int!**
    gpuErrchk(cudaMemcpyAsync( fpfb->d_filter_coeffs, fpfb->filter_coeffs, filter_size, cudaMemcpyHostToDevice ));

    // Allocate device memory for the other arrays
    if (flags & PFB_MALLOC_HOST_OUTPUT)
        gpuErrchk(cudaMallocHost( (void **)&(fpfb->vcs_data), fpfb->vcs_size ));

    if (flags & PFB_MALLOC_DEVICE_INPUT)
        gpuErrchk(cudaMalloc( (void **)&(fpfb->d_htr_data), fpfb->d_htr_size ));
    if (flags & PFB_MALLOC_DEVICE_OUTPUT)
        gpuErrchk(cudaMalloc( (void **)&(fpfb->d_vcs_data), fpfb->d_vcs_size ));
    gpuErrchk(cudaMalloc( (void **)&(fpfb->d_weighted_overlap_add), fpfb->weighted_overlap_add_size ));

    // Construct the cuFFT plan
    int rank     = 1;       // i.e. a 1D FFT
    int n        = fpfb->K; // size of each FFT
    int *inembed = NULL;    // Setting this to null makes all subsequent "data layout" parameters ignored
                            // to use the default layout (contiguous data in memory)

    fpfb->ninputs_per_cufft_batch = 32; // This seems to work, so I'll have to call cuFFT 256/32 = 8 times
    fpfb->cufft_batch_size = fpfb->nspectra_per_chunk * fpfb->ninputs_per_cufft_batch;

    cufftResult res = cufftPlanMany( &(fpfb->plan), rank, &n,
            inembed, 0, 0, NULL, 0, 0,  // <-- Here are all the ignored data layout parameters
            CUFFT_C2C, fpfb->cufft_batch_size );
    if (res != CUFFT_SUCCESS)
    {
        fprintf( stderr, "CUFFT error: Plan creation failed with error code %d\n", res );
        exit(EXIT_FAILURE);
    }
}

void free_forward_pfb( forward_pfb *fpfb )
/* Free host and device memory allocated in init_forward_pfb
 */
{
    cufftDestroy( fpfb->plan );
    gpuErrchk(cudaFreeHost( fpfb->filter_coeffs ));
    gpuErrchk(cudaFreeHost( fpfb->vcs_data ));
    gpuErrchk(cudaFreeHost( fpfb->i_output_idx ));
    gpuErrchk(cudaFree( fpfb->d_filter_coeffs ));
    gpuErrchk(cudaFree( fpfb->d_htr_data ));
    gpuErrchk(cudaFree( fpfb->d_vcs_data ));
    gpuErrchk(cudaFree( fpfb->d_i_output_idx ));
    gpuErrchk(cudaFree( fpfb->d_weighted_overlap_add ));
    free( fpfb );
}

void cu_forward_pfb( vcsbeam_context *vm, bool copy_result_to_host, logger *log )
/* The wrapper function that performs a forward PFB on the GPU
 */
{
    // Shorthand variables
    forward_pfb *fpfb = vm->fpfb;

    // Check that the input and output buffers have been set
    if (vm->v->buffer == NULL)
    {
        fprintf( stderr, "error: cu_forward_pfb: Input data buffer "
                "have not been set\n" );
        exit(EXIT_FAILURE);
    }

    if (copy_result_to_host && fpfb->vcs_data == NULL)
    {
        fprintf( stderr, "error: cu_forward_pfb: Output data buffer "
                "have not been set\n" );
        exit(EXIT_FAILURE);
    }

    for (fpfb->chunk = 0; fpfb->chunk < fpfb->chunks_per_second; fpfb->chunk++)
    {
        // Copy data to device
        logger_start_stopwatch( log, "upload", true );

        gpuErrchk(cudaMemcpy(
                    fpfb->d_htr_data,                                                 // to
                    (char *)vm->v->buffer + fpfb->chunk*fpfb->htr_stride,  // from
                    fpfb->d_htr_size,                                                 // how much
                    cudaMemcpyHostToDevice ));                                        // which direction

        logger_stop_stopwatch( log, "upload" );

        // Read lock can now be switched off
        if (fpfb->chunk == fpfb->chunks_per_second - 1)
            vm->v->locked = false;

        // PFB algorithm:
        // 1st step: weighted overlap add
        dim3 blocks( fpfb->nspectra_per_chunk, fpfb->I, fpfb->P );
        dim3 threads( fpfb->K );

        logger_start_stopwatch( log, "wola", true );

        // Set the d_weighted_overlap_add array to zeros
        gpuErrchk(cudaMemset( fpfb->d_weighted_overlap_add, 0, fpfb->weighted_overlap_add_size ));

        legacy_pfb_weighted_overlap_add<<<blocks, threads>>>( fpfb->d_htr_data, fpfb->d_filter_coeffs, fpfb->d_weighted_overlap_add );
        cudaDeviceSynchronize();
        gpuErrchk( cudaPeekAtLastError() );

        if (fpfb->flags & PFB_EMULATE_FPGA)
        {
            // Perform the weird rounding and demotion described in the appendix of McSweeney et al. (2020)
            dim3 blocks2( fpfb->nspectra_per_chunk, fpfb->I );
            dim3 threads2( fpfb->K );
            fpga_rounding_and_demotion<<<blocks2, threads2>>>( fpfb->d_weighted_overlap_add );
        }
        else
        {
            // Otherwise, just convert the ints to floats in preparation for the cuFFT
            dim3 blocks2( (2 * fpfb->nspectra_per_chunk * fpfb->I * fpfb->K) / 1024 ); // The factor of 2 is because each thread will deal with a single int, not a single int2
            dim3 threads2( 1024 );
            double scale = 1.0/16384.0; // equivalent to the ">> 14" operation applied in fpga_rounding_and_demotion()
            int2float<<<blocks2, threads2>>>( fpfb->d_weighted_overlap_add, scale );
        }
        cudaDeviceSynchronize();
        gpuErrchk( cudaPeekAtLastError() );

        logger_stop_stopwatch( log, "wola" );

        // 2nd step: FFT
        logger_start_stopwatch( log, "fft", true );

        int batch;
        for (batch = 0; batch < fpfb->I / fpfb->ninputs_per_cufft_batch; batch++)
        {
            cufftExecC2C(
                    fpfb->plan,
                    fpfb->d_weighted_overlap_add + batch*fpfb->cufft_batch_size*fpfb->K,
                    fpfb->d_weighted_overlap_add + batch*fpfb->cufft_batch_size*fpfb->K,
                    CUFFT_FORWARD );
        }
        cudaDeviceSynchronize();
        gpuErrchk( cudaPeekAtLastError() );

        logger_stop_stopwatch( log, "fft" );

        // 3rd step: packaging the result
        dim3 blocks3( fpfb->nspectra_per_chunk, fpfb->K );
        dim3 threads3( fpfb->I );

        logger_start_stopwatch( log, "pack", true );

        pack_into_recombined_format<<<blocks3, threads3>>>( fpfb->d_weighted_overlap_add,
                fpfb->d_vcs_data, fpfb->d_i_output_idx, fpfb->flags );
        cudaDeviceSynchronize();
        gpuErrchk( cudaPeekAtLastError() );

        logger_stop_stopwatch( log, "pack" );

        // Finally, copy the answer back to host memory, if requested
        if (copy_result_to_host)
        {
            logger_start_stopwatch( log, "download", true );

            gpuErrchk(cudaMemcpy(
                        (char *)fpfb->vcs_data + fpfb->chunk*fpfb->vcs_stride,   // to
                        fpfb->d_vcs_data,                                // from
                        fpfb->d_vcs_size,                                // how much
                        cudaMemcpyDeviceToHost ));                       // which direction

            logger_stop_stopwatch( log, "download" );
        }
    }
}

/**********************************
 * BACKWARDS (SYNTHESIS) FINE PFB *
 **********************************/

__global__ void ipfb_kernel(
    float *in_real, float *in_imag,
    float *ft_real, float *ft_imag,
    int ntaps, int npol, float *out )
/* This kernel computes the synthesis filter:

              1              K-1
   xhat[n] = --- SUM f[n-mM] SUM X_k[m] e^(2πjkn/K)
              K   m          k=0

   The sum over m is nominally over all integers, but in practice only
   involves a few terms because of the finiteness of the filter, f. To be
   precise, there are precisely ntaps non-zero values.

   X_k[m] represents the complex-valued inputs, in_real and in_imag.
   Every possible value of f[n]*e^(2πjkn/K) is provided in ft_real and
   ft_imag.

   K is the number of channels, and because this is a critically sampled
   PFB, M = K. We will also use P to mean the number of taps in the synthesis
   filter, and N = KP to mean the size of the filter.

   The polarisations are computed completely independently.

   And, of course, xhat[n] is represented by the out array.
 */
{
    // First, set a generic variable for this thread
    int idx = blockDim.x*blockIdx.x + threadIdx.x;

    // The polarisation for this thread is
    int pol = idx % npol;

    // and the time index of the output (i.e. the index for xhat) is
    int n = idx / npol;

    // Other constants we'll need are:
    int K = blockDim.x / npol;  // Total number of channels (should be 128)
    int M = K;
    int P = ntaps;
    int F = P*K;

    // Because we must have 0 <= n-mM < F, the smallest allowed value of m
    // is:
    int m0 = (n - F)/M + 1;

    // Initialise the output sample to zero
    float out_real = 0.0;
    float out_imag = 0.0;

    // Perform the double sum
    int m, k, f, tw, ft, i;
    for (m = m0; m < m0 + P; m++)
    {
        // With m now known, we can get the index for the filter
        f = n - m*M;

        //printf("n=%d, m=%d, f=%d\n", n, m, f);
        for (k = 0; k < K; k++)
        {
            // The index for the twiddle factor is
            tw = ((k+K/2)*n) % K;
            // (the extra K/2 identifies the middle channel as the DC bin)

            // The index into the ft (= filter/twiddle) array is
            ft = F*tw + f;

            // The "in" index (see cu_invert_pfb() for how the in[] arrays
            // were packed)
            // The fine channel time index, m, must be adjusted to ensure that
            // n=0 corresponds to the first full filter's worth of input samples
            i = npol*K*(m+P) + npol*k + pol;

            // Complex multiplication
            out_real += in_real[i] * ft_real[ft] -
                        in_imag[i] * ft_imag[ft];
            out_imag += in_real[i] * ft_imag[ft] +
                        in_imag[i] * ft_real[ft];
        }
    }

    // out[] includes both polarisations, at adjacent indices
    out[2*idx]   = out_real / K;
    out[2*idx+1] = out_imag / K;

    __syncthreads();
}

void cu_invert_pfb( cuDoubleComplex ****detected_beam, int file_no,
                        int npointing, int nsamples, int nchan, int npol,
                        int sizeof_buffer,
                        struct gpu_ipfb_arrays *g, float *data_buffer_vdif )
/* "Invert the PFB" by applying a resynthesis filter, using GPU
 * acceleration.
 *
 * This function expects "detected_beam" to be structured as follows:
 *
 *   detected_beam[2*nsamples][nchan][npol]
 *
 * Although detected_samples potentially contains 2 seconds' worth of data,
 * this function only inverts 1 second. The appropriate second is worked out
 * using file_no: if it is even, the first half of detected_beam is used,
 * if odd, the second half.
 *
 * The output of the inversion is packed back into data_buffer_vdif, a 1D
 * array whose ordering is as follows:
 *
 *   time, pol, complexity
 *
 * This ordering is suited for immediate output to the VDIF format.
 *
 * It is assumed that the inverse filter coefficients have already been loaded
 * to the GPU.
 */
{
    // Setup input values:
    // The starting sample index is "ntaps" places from the end of the second
    // half of detected_beam if the file number is even, and "ntaps" places
    // from the end of the first half of detected_beam if the file number is
    // odd.
    
    int start_s = (file_no % 2 == 0 ? 2*nsamples - g->ntaps : nsamples - g->ntaps);

    int p, s_in, s, ch, pol, i;
    for (p = 0; p < npointing; p++)
    for (s_in = 0; s_in < nsamples + g->ntaps; s_in++)
    {
        s = (start_s + s_in) % (2*nsamples);
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
            {
                // Calculate the index for in_real and in_imag;
                i = p    * npol * nchan * (nsamples + g->ntaps) +
                    s_in * npol * nchan +
                    ch   * npol +
                    pol;
                // Copy the data across - taking care of the file_no = 0 case
                // The s_in%(npol*nchan*nsamples) does this for each pointing
                if (file_no == 0 && (s_in%(npol*nchan*nsamples)) < g->ntaps)
                {
                    g->in_real[i] = 0.0;
                    g->in_imag[i] = 0.0;
                }
                else
                {
                    g->in_real[i] = cuCreal( detected_beam[p][s][ch][pol] );
                    g->in_imag[i] = cuCimag( detected_beam[p][s][ch][pol] );
                }
            }
        }
    }
    
    // Copy the data to the device
    gpuErrchk(cudaMemcpy( g->d_in_real, g->in_real, g->in_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( g->d_in_imag, g->in_imag, g->in_size, cudaMemcpyHostToDevice ));
    
    // Call the kernel
    if (npointing > 1)
    {
        fprintf( stderr, "error: PFB inversion currently only supports a single pointing\n" );
        exit(EXIT_FAILURE);
    }
    ipfb_kernel<<<nsamples, nchan*npol>>>( g->d_in_real, g->d_in_imag,
                                             g->d_ft_real, g->d_ft_imag,
                                             g->ntaps, npol, g->d_out );
    gpuErrchk( cudaPeekAtLastError() );
    cudaDeviceSynchronize();

    // Copy the result back into host memory
    gpuErrchk(cudaMemcpy( data_buffer_vdif, g->d_out, g->out_size, cudaMemcpyDeviceToHost ));
}


void cu_load_ipfb_filter( pfb_filter *filter, struct gpu_ipfb_arrays *g )
/* This function loads the inverse filter coefficients and the twiddle factors
   into GPU memory. If they were loaded separately (as floats), then the
   multiplication of the filter coefficients and the twiddle factors will be
   less precise than if a single array containing every combination of floats
   and twiddle factors is calculated in doubles, and then demoted to floats.
   Hence, this pre-calculation is done in this function before cudaMemcpy is
   called.

   The result is 2x 1D arrays loaded onto the GPU (one for real, one for imag)
   where the ith element is equal to

   result[i] = f[n] * exp(2πjk/K),
   n = i % N  (N is the filter size, "fil_size")
   k = i / N
   and K is the number of channels (nchan).

*/
{
    int ch, f, i;

    // Setup filter values:
    cuDoubleComplex ft; // pre-calculated filter coeffs times twiddle factor
    cuDoubleComplex cf; // temp variable for complex version of filter coeffs
    for (f = 0; f < filter->ncoeffs; f++)
    {
        cf = make_cuDoubleComplex( filter->coeffs[f], 0.0 );
        for (ch = 0; ch < filter->nchans; ch++)
        {
            i = filter->ncoeffs*ch + f;
            ft = cuCmul( filter->twiddles[ch], cf );
            g->ft_real[i] = cuCreal( ft );
            g->ft_imag[i] = cuCimag( ft );
        }
    }

    gpuErrchk(cudaMemcpy( g->d_ft_real, g->ft_real, g->ft_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( g->d_ft_imag, g->ft_imag, g->ft_size, cudaMemcpyHostToDevice ));
}


void malloc_ipfb( struct gpu_ipfb_arrays *g, pfb_filter *filter, int nsamples,
        int npol, int npointing )
{
    // Some shorthand variables:
    int ntaps = filter->ntaps;
    int nchan = filter->nchans;
    int fil_size = filter->ncoeffs;

    // Flatten the input array (detected_array) for GPU.
    // We only need one second's worth, plus 12 time samples tacked onto the
    // beginning (from the previous second)

    g->ntaps     = ntaps;
    g->in_size   = npointing * ((nsamples + ntaps) * nchan * npol) * sizeof(float);
    g->ft_size   = fil_size * nchan * sizeof(float);
    g->out_size  = npointing * nsamples * filter->nchans * npol * 2 * sizeof(float);

    // Allocate memory on the device
    gpuErrchk(cudaMalloc( (void **)&g->d_in_real, g->in_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_in_imag, g->in_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_ft_real, g->ft_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_ft_imag, g->ft_size ));

    gpuErrchk(cudaMalloc( (void **)&g->d_out, g->out_size ));

    // Allocate memory for host copies of the same
    g->in_real = (float *)malloc( g->in_size );
    g->in_imag = (float *)malloc( g->in_size );
    g->ft_real = (float *)malloc( g->ft_size );
    g->ft_imag = (float *)malloc( g->ft_size );

}


void free_ipfb( struct gpu_ipfb_arrays *g )
{
    // Free memory on host and device
    free( g->in_real );
    free( g->in_imag );
    free( g->ft_real );
    free( g->ft_imag );
    cudaFree( g->d_in_real );
    cudaFree( g->d_in_imag );
    cudaFree( g->d_ft_real );
    cudaFree( g->d_ft_imag );
    cudaFree( g->d_out );
}
