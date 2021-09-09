/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>

extern "C" {
#include "pfb.h"
#include "filter.h"
#include "jones.h"
}

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
        int *filter_coeffs, cuFloatComplex *weighted_overlap_array )
{
    // Parse the block and thread idxs:
    //   <<<(nspectra,I),(K,P)>>>
    // where nspectra is the number of output spectra,
    //       I        is the number of RF inputs,
    //       K        is the size of the output spectrum
    //       P        is the number of taps
    // ...and put everything into the mathematical notation used in McSweeney
    // et al. (2020) (see equation in comments above)
    int              m        = blockIdx.x;
    int              nspectra = gridDim.x;
    int              I        = gridDim.y;
    int              i        = blockIdx.y;
    int              K        = blockDim.x;
    int              M        = K; // This enforces a critical sampled PFB
    //int              P        = blockDim.y;
    int              n        = threadIdx.x;
    int              p        = threadIdx.y;

    int             *h = filter_coeffs;
    char2           *x = indata;
    cuFloatComplex  *b = weighted_overlap_array;

    // Use shared memory as a temporary workspace for preparing the b array
    // For one block, this should have K elements
    extern __shared__ int2 bint[];

    // Let the first tap (p=0) have the responsibility of initialising the
    // bint array to zeros
    if (p == 0)
    {
        bint[n].x = 0;
        bint[n].y = 0;
    }
    __syncthreads();

    // Now calculate the index into the various arrays that also
    // takes into account the fact that these arrays contain all RF inputs.
    // MEMO TO SELF: My current going theory is that I don't have to do any
    // re-ordering of the antennas, as that is dealt with elsewhere.
    int h_idx = K*p - n;
    int x_idx = vMWAX_IDX(n + m*M - K*p, i, I);
    int b_idx = (K*i + m)*nspectra + n; // This puts each set of K samples to
                                        // be FFT'd in a contiguous memory block

    // Now perform the weighted overlap-add operation
    int   hval = h[h_idx];
    char2 xval = x[x_idx];

    atomicAdd( &bint[n].x, hval*(int)xval.x );
    atomicAdd( &bint[n].y, hval*(int)xval.y );

    __syncthreads();

    // In keeping with the original FPGA implementation, the result now needs to
    // be demoted and rounded. Only one tap needs to do this
    int X, Y; // To avoid too many shared memory accesses
    if (p == 0)
    {
        X = bint[n].x;
        Y = bint[n].y;

        // Rounding:
        if (X > 0)  X += 0x2000;
        if (Y > 0)  Y += 0x2000;

        // Demotion:
        X >>= 14;
        Y >>= 14;

        // Promote the result to doubles and put it in the b array in global memory
        // in preparation for being FFTed
        b[b_idx] = make_cuFloatComplex( (double)X, (double)Y );
    }

    __syncthreads();
}

__global__ void pack_into_recombined_format( cuFloatComplex *ffted, uint8_t *outdata )
/* This is the final step in the forward fine PFB algorithm that emulates what
   was implemented on the MWA FPGAs in Phase 1 & 2 (see above for details).
   At this point, the FFTED array contains the Fourier-transformed data that
   already represents the final channelisation. All that remains to be done is
   to pack it into the same format as the VCS recombined data.

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
    int nspectra = gridDim.x;
    int K        = gridDim.y;
    int m        = blockIdx.x;
    int k        = blockIdx.y;
    int i        = threadIdx.x;
    int I        = blockIdx.x;

    cuFloatComplex  *b = ffted;
    uint8_t         *X = outdata;

    // Calculate the idxs into b and X
    int b_idx = (K*i + m)*nspectra + k;
    int X_idx = v_IDX(m, k, i, K, I);

    // Pull the values to be manipulated into register memory (because the
    // packing macro below involves a lot of repetition of the arguments)
    double re = b[b_idx].x;
    double im = b[b_idx].y;

    // Put the packed value back into global memory at the appropriate place
    X[X_idx] = PACK_NIBBLES(re, im);

    __syncthreads();
}

forward_pfb *init_forward_pfb( MetafitsMetadata *obs_metadata,
        char2 *htr_data, char2 *htr_data_extended, uint8_t *vcs_data,
        pfb_filter *filter )
/* Create and initialise a forward_pfb struct.

   Inputs:
     OBS_METADATA      - mwalib metadata struct
     HTR_DATA          - pointer to host memory to be PFB'd
     HTR_DATA_EXTENDED - pointer to extended host memory to be PFB'd
                         (will be tacked onto the end of

   Output:
     FPFB         - Pointer to struct to be initialised
 */
{
    // Create the struct in memory
    forward_pfb *fpfb = (forward_pfb *)malloc( sizeof(forward_pfb) );

    // Host memory is assumed to be allocated
    fpfb->htr_data          = htr_data;
    fpfb->htr_data_extended = htr_data_extended;
    fpfb->vcs_data          = vcs_data;
}

void free_forward_pfb( forward_pfb *pfb )
// Free the memory allocated in init_forward_pfb
{
    free( fpfb );
}

void cu_forward_pfb_fpga_version( forward_pfb *fpfb, bool copy_result_to_host )
/* The wrapper function that performs the forward PFB algorithm as originally
   implemented on the FPGAs for MWA Phases 1 & 2.
   A cuFFT plan must already have been made, via make_forward_pfb_fpga_fft_plan().
 */
{
    // Copy data to device
    gpuErrchk(cudaMemcpy( fpfb->d_htr_data, fpfb->htr_data, fpfb->htr_size, cudaMemcpyHostToDevice ));
    if (fpfb->htr_data_extended != NULL)
        gpuErrchk(cudaMemcpy( fpfb->d_htr_data + fpfb->htr_size, fpfb->htr_data_extended, fpfb->htr_extended_size, cudaMemcpyHostToDevice ));

    // PFB algorithm:
    // 1st step: weighted overlap add
    dim3 blocks( fpfb->nspectra, fpfb->I );
    dim3 threads( fpfb->K, fpfb->P );

    legacy_pfb_weighted_overlap_add<<<blocks, threads>>>( fpfb->d_htr_data, fpfb->d_filter_coeffs, fpfb->d_weighted_overlap_add );
    gpuErrchk( cudaPeekAtLastError() );

    // 2nd step: FFT
    cufftExecC2C( fpfb->plan, fpfb->d_weighted_overlap_add, fpfb->d_weighted_overlap_add, CUFFT_FORWARD );

    // 3rd step: packaging the result
    dim3 blocks2( fpfb->nspectra, fpfb->K );
    dim3 threads2( fpfb->I );

    pack_into_recombined_format<<<blocks2, threads2>>>( fpfb->d_weighted_overlap_add, fpfb->d_vcs_data );
    gpuErrchk( cudaPeekAtLastError() );

    // Finally, copy the answer back to host memory, if requested
    if (copy_result_to_host)
    {
        gpuErrchk(cudaMemcpy( fpfb->vcs_data, fpfb->d_vcs_data, fpfb->vcs_size, cudaMemcpyDeviceToHost ));
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
    for (f = 0; f < filter->size; f++)
    {
        cf = make_cuDoubleComplex( filter->coeffs[f], 0.0 );
        for (ch = 0; ch < filter->nchans; ch++)
        {
            i = filter->size*ch + f;
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
    int fil_size = filter->size;

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
