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
#include "gpu_fft.hpp"
#include "gpu_macros.h"

#include <mwalib.h>

#include "vcsbeam.h"

inline void gpuAssert(gpuError_t code, const char *file, int line, bool abort=true)
{
    /* Wrapper function for GPU/CUDA error handling. Every CUDA call goes through
       this function. It will return a message giving your the error string,
       file name and line of the error. Aborts on error. */

    if (code != 0)
    {
        fprintf(stderr, "GPUAssert:: %s - %s (%d)\n", gpuGetErrorString(code), file, line);
        if (abort)
        {
            exit(code);
        }
    }
}

// define a macro for accessing gpuAssert
// #define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__, true);}


/**
 * \file pfb.cu
 *
 * [McSweeney2020]: https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/mwa-tiedarray-processing-iii-microsecond-time-resolution-via-a-polyphase-synthesis-filter/C76837FE84A1DCFAA696C9AC10C44F2D
 * [MWAXHTR]: https://wiki.mwatelescope.org/display/MP/MWA+High+Time+Resolution+Voltage+Capture+System
 *
 * # Forward (analysis) fine PFB
 *
 * The following kernels are part of an "offline" version of the "fine PFB"
 * algorithm that was implemented on the FPGAs of Phase 1 & 2 of the MWA. As
 * described in [McSweeney et al. (2020)][McSweeney2020], this algorithm is a
 * version of the "weighted overlap-add" algorithm (see their Eq. (3)):
 * \f[
 *     X_k[m] \sum_n^{K-1} b_m[n] e^{-2\pi jkn/K},
 * \f]
 * where
 * \f[
 *     b_m[n] = \sum_\rho^{P-1} h[K\rho - n] x[n + mM - K\rho].
 * \f]
 *
 * A full description of these symbols is given in the reference, but here
 * is a summary:
 *
 * | Symbol | Meaning | Variable |
 * | :----: | :------ | :------- |
 * | \f$x\f$ | The input data                  | `indata`                 |
 * | \f$X\f$ | The output data                 | `outdata`                |
 * | \f$h\f$ | The filter coefficients         | `filter_coeffs`          |
 * | \f$P\f$ | The number of taps              | `ntaps`                  |
 * | \f$b\f$ | The weighted overlap-add array  | `weighted_overlap_array` |
 * | \f$\text{nspectra}\f$ | The number of output spectra |               |
 * | \f$I\f$ | The number of RF inputs         |                          |
 * | \f$K\f$ | The size of the output spectrum |                          |
 * | \f$P\f$ | The number of taps              |                          |
 *
 * Note that this nomenclature differs from that used elsewhere in this
 * documentation.
 *
 * The algorithm is broken up into three parts, the third of which is
 * optional:
 *  1. the weighted overlap-add, which forms `b` from `x` and `h`,
 *  2. the FFT, which is implemented using the cuFFT library, and
 *  3. the demotion and packaging of the result into the
 *     output format, which is described in the Appendix of
 *     [McSweeney et al. (2020)][McSweeney2020]
 *
 * ## Notes:
 *
 *  - `indata` is expected to have the same data layout as delivered by mwalib's
 *    mwalib_voltage_context_read_second() function when applied to MWAX high
 *    time resolution data, whose format description can be found
 *    [here][MWAXHTR]. The samples are organised according to the `vMWAX_IDX`
 *    macro.
 *  - `outdata` will have the same data layout as "recombined" legacy VCS data,
 *    as per the `v_IDX` macro.
 *  - Each thread will operate on one RF input (i.e. antenna/pol combination)
 *    and generate the spectrum for a single "fine-channelised" time step
 *    (the `m` index). The `weighted overlap-add` array (`b`) also resides in
 *    device memory, in order that it can make use of cuFFT.
 */

/**
 * Performs the weighted overlap-add part of the PFB algorithm.
 *
 * @param[in] indata The input data, \f$x[n]\f$,
 *                   with layout equivalent to the buffer populated by
 *                   mwalib (see [the MWAX voltage format][MWAXHTR])
 * @param[in] filter_coeffs The PFB filter coefficients, \f$h[n]\f$
 * @param[out] weighted_overlap_array The result of the weighted
 *             overlap-add operation, \f$b[n]\f$
 *
 * The weighted overlap-add part of the PFB algorithm is the equation
 * \f[
 *     b_m[n] = \sum_\rho^{P-1} h[K\rho - n] x[n + mM - K\rho].
 * \f]
 * See above for the meaning of the terms in this equation.
 * 
 * The expected thread configuration is
 * \f$ \langle\langle\langle (\text{nspectra},I,P), K \rangle\rangle\rangle \f$
 *
 * \todo Add data layout information to this docstring.
 */
__global__ void vmWOLA_kernel( char2 *indata,
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

/**
 * CUDA kernel for performing the rounding and demotion step of the forward
 * PFB.
 *
 * @param[in,out] data The data to be rounded and demoted.
 *
 * The rounding and demotion performed here is designed to emulate that
 * done in the FPGAs of the legacy (Phases 1 & 2) MWA system.
 * Specifically, it does the following:
 * \f[
 *     \begin{cases}
 *         \left\lfloor \frac{x}{2^{14}} + 0.5 \right\rfloor, & x > 0, \\
 *         \left\lfloor \frac{x}{2^{14}}       \right\rfloor, & x \le 0.
 *     \end{cases}
 * \f]
 The final numbers are represented as 32-bit floats.
 *
 * The expected thread configuration is
 * \f$\langle\langle\langle (\text{nspectra} \times I), K
 * \rangle\rangle\rangle\f$
 */
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
    void* ptr = &(b[b_idx].x);
    float *fx = (float *)ptr;
    
    ptr = &(b[b_idx].y);
    float *fy = (float *)ptr;

    *fx = (float)X;
    *fy = (float)Y;
}

/**
 * Convert 32-bit integers to 32-bit floats.
 *
 * @param data The data to be converted
 * @param scale A scaling factor to apply (set to 1 for no effect)
 *
 * Each thread handles a single integer, so the thread configuration can be
 * chosen to suit the size of `data`.
 *
 * This function acts as an alternative to fpga_rounding_and_demotion(),
 * in order to preserve the full precision of the PFB output.
 */
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

/**
 * CUDA kernel for packing the PFB'ed data into legacy recombined
 * layout/format.
 *
 * @param[in] ffted The already-FFT'ed PFB data, ready for packing
 * @param[out] outdata The output data buffer
 * @param i_idx An array of indexes for the desired ordering of the RF inputs
 *              (i.e. antenna/polarisation combinations)
 * @param flags PFB configuration settings
 *
 * This is the final step in the forward fine PFB algorithm.
 * At this point, the FFTED array contains the Fourier-transformed data that
 * already represents the final channelisation. All that remains to be done is
 * to pack it into the same layout (and optionally also the same format) as
 * the VCS recombined data.
 *
 * The relevant values for `flags` are (see vmInitForwardPFB() for the full
 * table):
 *
 * | Flag name | Description |
 * | :-------- | :---------- |
 * | `PFB_COMPLEX_INT4`    | Typecast the output into (4+4)-bit complex integers |
 * | `PFB_COMPLEX_FLOAT64` | Typecast the output into (64+64)-bit complex floats |
 *
 * The expected thread configuration is
 * \f$\langle\langle\langle(\text{nspectra},K),I\rangle\rangle\rangle\f$.
 */
__global__ void pack_into_recombined_format( gpuFloatComplex *ffted, void *outdata, int *i_idx, pfb_flags flags )
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

    gpuFloatComplex  *b = ffted;

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
            X[X_idx] = PACK_NIBBLES(re, im);
        else
            X[X_idx] = PACK_NIBBLES(im, re);
    }
    else // Currently, default is gpuDoubleComplex
    {
        gpuDoubleComplex *X = (gpuDoubleComplex *)outdata;
        if (flags & PFB_IMAG_PART_FIRST)
            X[X_idx] = make_gpuDoubleComplex( re, im );
        else
            X[X_idx] = make_gpuDoubleComplex( im, re );
    }

    __syncthreads();
}

/**
 * Create and initialise a forward_pfb struct.
 *
 * @param vm The VCSBeam context struct
 * @param M The "stride" of the PFB
 * @param flags PFB configuration settings
 *
 * Setting \f$M = K\f$, where \f$K\f$ is the value of
 * `vm&rarr;analysis_filter&rarr;nchans`, will make this a
 * "critically sampled PFB"; setting \f$M > K\f$, an "oversampled PFB".
 *
 * A filter must already be loaded into `vm&rarr;analysis_filter`.
 * **WARNING! The filter will be forcibly typecast to int!!**
 *
 * The GPU processing will be divided up into `vm&rarr;chunks_per_second`,
 * chunks per second, but if this number does not divide the number
 * of samples (10000) evenly, then it is incremented until it does.
 *
 * The available flags are given in the following table:
 * | Flag name | Description |
 * | :-------- | :---------- |
 * | `PFB_MALLOC_HOST_INPUT`    | Allocate memory on CPU for the input        |
 * | `PFB_MALLOC_HOST_OUTPUT`   | Allocate memory on CPU for the output       |
 * | `PFB_MALLOC_DEVICE_INPUT`  | Allocate memory on GPU for the input        |
 * | `PFB_MALLOC_DEVICE_OUTPUT` | Allocate memory on GPU for the output       |
 * | `PFB_MALLOC_ALL`           | Synonym for <code>PFB_MALLOC_HOST_INPUT \| PFB_MALLOC_HOST_OUTPUT \| PFB_MALLOC_DEVICE_INPUT \| PFB_MALLOC_DEVICE_OUTPUT</code> |
 * | `PFB_COMPLEX_INT4`         | Typecast the output into (4+4)-bit complex integers |
 * | `PFB_COMPLEX_FLOAT64`      | Typecast the output into (64+64)-bit complex floats |
 * | `PFB_IMAG_PART_FIRST`      | Place the imaginary part of the output in the first (i.e. most significant) position |
 * | `PFB_EMULATE_FPGA`         | Perform the (asymmetric) rounding and demotion step in exactly the same way as the original (Phase 1 & 2) MWA FPGAs |
 * | `PFB_SMART`                | Synonym for <code>PFB_MALLOC_HOST_INPUT \| PFB_MALLOC_DEVICE_INPUT \| PFB_MALLOC_DEVICE_OUTPUT \| PFB_EMULATE_FPGA \| PFB_COMPLEX_INT4</code> |
 * | `PFB_FULL_PRECISION`       | Synonym for <code>PFB_MALLOC_HOST_INPUT \| PFB_MALLOC_DEVICE_INPUT \| PFB_MALLOC_DEVICE_OUTPUT \| PFB_COMPLEX_FLOAT64</code> |
 */
void vmInitForwardPFB( vcsbeam_context *vm, int M, pfb_flags flags )
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
    (gpuMallocHost( (void **)&(fpfb->i_output_idx),   fpfb->I*sizeof(int) ));
    (gpuMalloc(     (void **)&(fpfb->d_i_output_idx), fpfb->I*sizeof(int) ));

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

    (gpuMemcpyAsync( fpfb->d_i_output_idx, fpfb->i_output_idx, fpfb->I*sizeof(int), gpuMemcpyHostToDevice ));

    // Work out the sizes of the various arrays
    while (fpfb->nspectra % vm->chunks_per_second != 0)
        vm->chunks_per_second++;
    fpfb->nspectra_per_chunk = fpfb->nspectra / vm->chunks_per_second;

    int num_voltage_blocks_per_chunk = vm->vcs_metadata->num_voltage_blocks_per_second / vm->chunks_per_second;
    fpfb->d_htr_size = (num_voltage_blocks_per_chunk + 1) * vm->vcs_metadata->voltage_block_size_bytes;

    if (flags & PFB_COMPLEX_INT4)
    {
        fpfb->vcs_size = fpfb->nspectra * vm->obs_metadata->num_rf_inputs * fpfb->K * sizeof(uint8_t);
        vm->datatype = VM_INT4;
    }
    else // i.e., default is gpuDoubleComplex
    {
        fpfb->vcs_size = fpfb->nspectra * vm->obs_metadata->num_rf_inputs * fpfb->K * sizeof(gpuDoubleComplex);
        vm->datatype = VM_DBL;
    }

    fpfb->d_vcs_size = fpfb->vcs_size / vm->chunks_per_second;

    fpfb->htr_stride = (vm->vcs_metadata->num_voltage_blocks_per_second*vm->vcs_metadata->voltage_block_size_bytes)/vm->chunks_per_second;
    fpfb->vcs_stride = fpfb->d_vcs_size;

    fpfb->char2s_per_second =
        (vm->vcs_metadata->num_voltage_blocks_per_second *
         vm->vcs_metadata->voltage_block_size_bytes) / sizeof(char2);
    fpfb->bytes_per_block = vm->vcs_metadata->voltage_block_size_bytes;

    fpfb->weighted_overlap_add_size = vm->v->buffer_size * (sizeof(gpuFloatComplex) / sizeof(char2)) / vm->chunks_per_second;
    size_t filter_size = vm->analysis_filter->ncoeffs * sizeof(int);

    // Allocate memory for filter and copy across the filter coefficients,
    // casting to int
    (gpuMallocHost( (void **)&(fpfb->filter_coeffs),   filter_size ));
    (gpuMalloc(     (void **)&(fpfb->d_filter_coeffs), filter_size ));
    for (i = 0; i < vm->analysis_filter->ncoeffs; i++)
        fpfb->filter_coeffs[i] = (int)vm->analysis_filter->coeffs[i]; // **WARNING! Forcible typecast to int!**
    (gpuMemcpyAsync( fpfb->d_filter_coeffs, fpfb->filter_coeffs, filter_size, gpuMemcpyHostToDevice ));

    // Allocate device memory for the other arrays
    if (flags & PFB_MALLOC_HOST_OUTPUT)
    {
        gpuMallocHost( (void **)&(fpfb->vcs_data), fpfb->vcs_size );
    }
    else
        fpfb->vcs_data = NULL;

    if (flags & PFB_MALLOC_DEVICE_INPUT)
        (gpuMalloc( (void **)&(fpfb->d_htr_data), fpfb->d_htr_size ));
    if (flags & PFB_MALLOC_DEVICE_OUTPUT)
        (gpuMalloc( (void **)&(fpfb->d_vcs_data), fpfb->d_vcs_size ));
    (gpuMalloc( (void **)&(fpfb->d_weighted_overlap_add), fpfb->weighted_overlap_add_size ));

    // Construct the cuFFT plan
    int rank     = 1;       // i.e. a 1D FFT
    int n        = fpfb->K; // size of each FFT
    int *inembed = NULL;    // Setting this to null makes all subsequent "data layout" parameters ignored
                            // to use the default layout (contiguous data in memory)

    fpfb->ninputs_per_cufft_batch = 32; // This seems to work, so I'll have to call cuFFT 256/32 = 8 times
    fpfb->cufft_batch_size = fpfb->nspectra_per_chunk * fpfb->ninputs_per_cufft_batch;

    gpufftResult res = gpufftPlanMany( &(fpfb->plan), rank, &n,
            inembed, 0, 0, NULL, 0, 0,  // <-- Here are all the ignored data layout parameters
            GPUFFT_C2C, fpfb->cufft_batch_size );
    if (res != GPUFFT_SUCCESS)
    {
        fprintf( stderr, "GPUFFT error: Plan creation failed with error code %d\n", res );
        exit(EXIT_FAILURE);
    }

    vm->fine_sample_rate = fpfb->nspectra;
    vm->nfine_chan       = fpfb->K;
}

/**
 * Frees host and device memory allocated in vmInitForwardPFB().
 *
 * @param fpfb The `forward_pfb` object to be freed.
 *
 * After freeing the member variables of `fpfb`, the struct itself is freed.
 */
void vmFreeForwardPFB( forward_pfb *fpfb )
{
    gpufftDestroy( fpfb->plan );
    (gpuHostFree( fpfb->filter_coeffs ));
    (gpuHostFree( fpfb->vcs_data ));
    (gpuHostFree( fpfb->i_output_idx ));
    (gpuFree( fpfb->d_filter_coeffs ));
    (gpuFree( fpfb->d_htr_data ));
    (gpuFree( fpfb->d_vcs_data ));
    (gpuFree( fpfb->d_i_output_idx ));
    (gpuFree( fpfb->d_weighted_overlap_add ));
    free( fpfb );
}

/**
 * Copies a chunk of PFB data from CPU to GPU.
 *
 * A diagram of the memory layout of the copy operation:
 * ```
 *          <--hs->
 * Host:   |-------|-------|-------| ... |-------| ...
 * chunk:  0       1       2       3     c
 *                                       ^ from
 *
 *          <--ds-->
 * Device: |--------|
 *         ^ to
 * ```
 * where
 *  - **from** = `vm&rarr;v&rarr;buffer + c*hs`
 *  - **to** = `vm&rarr;fpfb&rarr;d_htr_data`
 *  - **ds** = `vm&rarr;fpfb&rarr;d_htr_size`
 *  - **hs** = `vm&rarr;fpfb&rarr;htr_stride`
 *  - **c** = `vm&rarr;chunk_to_load % vm&rarr;chunks_per_second`
 *
 * Note that **ds** is not necessarily equal to **hs**, and that **ds** bytes are
 * copied.
 */
void vmUploadForwardPFBChunk( vcsbeam_context *vm )
{
    logger_start_stopwatch( vm->log, "upload", false );

    int chunk = vm->chunk_to_load % vm->chunks_per_second;

    gpuMemcpy(
            vm->fpfb->d_htr_data,                                // to
            (char *)vm->v->buffer + chunk*vm->fpfb->htr_stride,  // from
            vm->fpfb->d_htr_size,                                // how much
            gpuMemcpyHostToDevice );                            // which direction

    logger_stop_stopwatch( vm->log, "upload" );

    // If it's the last chunk of the second, read lock can be switched off
    if (chunk == vm->chunks_per_second - 1)
        vm->v->locked = false;
}

/**
 * Calls the vmWOLA_kernel() kernel for the forward PFB data.
 *
 * @todo Keep this as a "pure" forward_pfb function, and make a separate
 *       "vm" function that's exposed.
 */
void vmWOLAChunk( vcsbeam_context *vm )
{
    // Shorthand variable
    forward_pfb *fpfb = vm->fpfb;

    dim3 blocks( fpfb->nspectra_per_chunk, fpfb->I, fpfb->P );
    dim3 threads( fpfb->K );

    logger_start_stopwatch( vm->log, "pfb-wola", false );

    // Set the d_weighted_overlap_add array to zeros
    (gpuMemset( fpfb->d_weighted_overlap_add, 0, fpfb->weighted_overlap_add_size ));

    vmWOLA_kernel<<<blocks, threads>>>( fpfb->d_htr_data, fpfb->d_filter_coeffs, fpfb->d_weighted_overlap_add );
    gpuDeviceSynchronize();
    ( gpuPeekAtLastError() );

    logger_stop_stopwatch( vm->log, "pfb-wola" );
}

/**
 * Calls the fpga_rounding_and_demotion() kernel for the forward PFB data.
 *
 * @todo Keep this as a "pure" forward_pfb function, and make a separate
 *       "vm" function that's exposed.
 */
void vmFPGARoundingChunk( vcsbeam_context *vm )
{
    // Shorthand variable
    forward_pfb *fpfb = vm->fpfb;

    logger_start_stopwatch( vm->log, "pfb-round", false );

    if (fpfb->flags & PFB_EMULATE_FPGA)
    {
        // Perform the weird rounding and demotion described in the appendix of McSweeney et al. (2020)
        dim3 blocks( fpfb->nspectra_per_chunk, fpfb->I );
        dim3 threads( fpfb->K );
        fpga_rounding_and_demotion<<<blocks, threads>>>( fpfb->d_weighted_overlap_add );
    }
    else
    {
        // Otherwise, just convert the ints to floats in preparation for the gpuFFT
        dim3 blocks( (2 * fpfb->nspectra_per_chunk * fpfb->I * fpfb->K) / 1024 ); // The factor of 2 is because each thread will deal with a single int, not a single int2
        dim3 threads( 1024 );
        double scale = 1.0/16384.0; // equivalent to the ">> 14" operation applied in fpga_rounding_and_demotion()
        int2float<<<blocks, threads>>>( fpfb->d_weighted_overlap_add, scale );
    }
    gpuDeviceSynchronize();
    ( gpuPeekAtLastError() );

    logger_stop_stopwatch( vm->log, "pfb-round" );
}

/**
 * Executes the (CUDA) FFT on the forward PFB data.
 *
 * @todo Keep this as a "pure" forward_pfb function, and make a separate
 *       "vm" function that's exposed.
 */
void vmFFTChunk( vcsbeam_context *vm )
{
    // Shorthand variable
    forward_pfb *fpfb = vm->fpfb;

    logger_start_stopwatch( vm->log, "pfb-fft", false );

    int batch;
    for (batch = 0; batch < fpfb->I / fpfb->ninputs_per_cufft_batch; batch++)
    {
        gpufftExecC2C(
                fpfb->plan,
                fpfb->d_weighted_overlap_add + batch*fpfb->cufft_batch_size*fpfb->K,
                fpfb->d_weighted_overlap_add + batch*fpfb->cufft_batch_size*fpfb->K,
                GPUFFT_FORWARD );
    }
    gpuDeviceSynchronize();
    ( gpuPeekAtLastError() );

    logger_stop_stopwatch( vm->log, "pfb-fft" );
}

/**
 * Calls the pack_into_recombined_format() kernel for the forward PFB data.
 *
 * @todo Keep this as a "pure" forward_pfb function, and make a separate
 *       "vm" function that's exposed.
 */
void vmPackChunk( vcsbeam_context *vm )
{
    // Shorthand variable
    forward_pfb *fpfb = vm->fpfb;

    dim3 blocks( fpfb->nspectra_per_chunk, fpfb->K );
    dim3 threads( fpfb->I );

    logger_start_stopwatch( vm->log, "pfb-pack", false );

    pack_into_recombined_format<<<blocks, threads>>>( fpfb->d_weighted_overlap_add,
            fpfb->d_vcs_data, fpfb->d_i_output_idx, fpfb->flags );
    gpuDeviceSynchronize();
    ( gpuPeekAtLastError() );

    logger_stop_stopwatch( vm->log, "pfb-pack" );
}

/**
 * Copies a chunk of PFB data from GPU to CPU.
 *
 * A diagram of the memory layout of the copy operation:
 * ```
 *          <--ds-->
 * Device: |--------|
 *         ^ from
 *
 *          <--hs->
 * Host:   |-------|-------|-------| ... |-------| ...
 * chunk:  0       1       2       3     c
 *                                       ^ to
 * ```
 * where
 *  - **from** = `vm&rarr;fpfb&rarr;d_vcs_data`
 *  - **to** = `vm&rarr;fpfb&rarr;vcs_data + c*hs`
 *  - **ds** = `vm&rarr;fpfb&rarr;d_vcs_size`
 *  - **hs** = `vm&rarr;fpfb&rarr;vcs_stride`
 *  - **c** = `vm&rarr;chunk_to_load % vm&rarr;chunks_per_second`
 *
 * Note that **ds** is not necessarily equal to **hs**, and that **ds** bytes are
 * copied.
 */
void vmDownloadForwardPFBChunk( vcsbeam_context *vm )
{
    // Shorthand variable
    forward_pfb *fpfb = vm->fpfb;

    logger_start_stopwatch( vm->log, "download", false );

    int chunk = vm->chunk_to_load % vm->chunks_per_second;

    (gpuMemcpy(
                (char *)fpfb->vcs_data + chunk*fpfb->vcs_stride,   // to
                fpfb->d_vcs_data,                                // from
                fpfb->d_vcs_size,                                // how much
                gpuMemcpyDeviceToHost ));                       // which direction

    logger_stop_stopwatch( vm->log, "download" );
}

/**
 * The wrapper function that performs a forward PFB on the GPU.
 *
 * For each chunk of data, it calls, in order,
 *   1. vmUploadForwardPFBChunk()
 *   2. vmWOLAChunk()
 *   3. vmFPGARoundingChunk()
 *   4. vmFFTChunk()
 *   5. vmPackChunk()
 *
 * If a buffer (`vm&rarr;fpfb&rarr;vcs_data`) has been prepared, then
 * vmDownloadForwardPFBChunk() is called.
 */
void vmExecuteForwardPFB( vcsbeam_context *vm )
{
    // Shorthand variables
    forward_pfb *fpfb = vm->fpfb;

    // Check that the input buffer has been set
    if (vm->v->buffer == NULL)
    {
        fprintf( stderr, "error: vmExecuteForwardPFB: Input data buffer "
                "have not been set\n" );
        exit(EXIT_FAILURE);
    }

    logger_start_stopwatch( vm->log, "pfb", true );

    int chunk;
    for (chunk = 0; chunk < vm->chunks_per_second; chunk++)
    {
        // Copy data to device
        vmUploadForwardPFBChunk( vm );

        // PFB algorithm:
        // 1) weighted overlap add
        vmWOLAChunk( vm );

        // 2) Rounding/demotion/scaling
        vmFPGARoundingChunk( vm );

        // 3) FFT
        vmFFTChunk( vm );

        // 4) Packaging the result
        vmPackChunk( vm );

        // Finally, copy the answer back to host memory, if set up to do so
        if (fpfb->vcs_data != NULL)
            vmDownloadForwardPFBChunk( vm );

        // Increment chunk counter
        vm->chunk_to_load++;
    }

    logger_stop_stopwatch( vm->log, "pfb" );
}

/**
 * Writes the result of the forward PFB operation to file.
 *
 * This function writes the contents of `vm&rarr;fpfb&rarr;vcs_data` to a file
 * whose name is returned by the function vmGetLegacyVoltFilename().
 */
void vmWritePFBOutputToFile( vcsbeam_context *vm )
{
    logger_start_stopwatch( vm->log, "write", true );

    char filename[128];
    vmGetLegacyVoltFilename( vm,
            vm->coarse_chan_idxs_to_process[0],
            vm->gps_seconds_to_process[vm->current_gps_idx-1],
            filename );

    FILE *f = fopen( filename, "w" );
    fwrite( vm->fpfb->vcs_data, vm->fpfb->vcs_size, sizeof(uint8_t), f );
    fclose( f );

    logger_stop_stopwatch( vm->log, "write" );
}

/**
 * \file pfb.cu
 *
 * # Backwards (synthesis) fine PFB
 *
 * The backwards, inverse, or synthesis PFB is implemented in a single CUDA
 * kernel, so see ipfb_kernel() for more information.
 */

/**
 * CUDA kernel implementing the synthesis PFB.
 *
 * @param[in]  in_real ...
 * @param[in]  in_imag ...
 * @param[in]  ft_real ...
 * @param[in]  ft_imag ...
 * @param      ntaps ...
 * @param      npol ...
 * @param[out] out ...
 *
 * The backwards/inverse/synthesis filter is
 * \f[
 *     \hat{x}[n] = \frac{1}{K} \sum_m f[n - mM]
 *                  \sum_{k=0}^{K-1} X_k[m]\,e^{2\pi jkn/K}
 * \f]
 *
 * The sum over `m` is nominally over all integers, but in practice only
 * involves a few terms because of the finiteness of the filter, `f`. To be
 * precise, there are precisely `ntaps` non-zero values.
 *
 * \f$X_k[m]\f$ represents the complex-valued inputs, `in_real` and `in_imag`.
 * Every possible value of \f$f[n]\,e^{2\pi jkn/K}\f$ is provided in `ft_real`
 * and `ft_imag`.
 *
 * `K` is the number of channels, and because this is a critically sampled
 * PFB, \f$M = K\f$. We will also use `P` to mean the number of taps in the
 * synthesis filter, and \f$N = KP\f$ to mean the size of the filter.
 *
 * The polarisations are computed completely independently.
 *
 * \f$\hat{x}[n]\f$ is represented by the `out` array.
 *
 * \todo Revamp this function to use cuComplex.
 */
__global__ void ipfb_kernel(
    float *in_real, float *in_imag,
    float *ft_real, float *ft_imag,
    int ntaps, int npol, float *out )
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

    // Because we must have 0 <= n-mM < F, the smallest allowed value of m is:
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
            i = (m+P) * npol * K + \
                k     * npol     + \
                pol;

            // Complex multiplication
            out_real += in_real[i] * ft_real[ft] -
                        in_imag[i] * ft_imag[ft];
            out_imag += in_real[i] * ft_imag[ft] +
                        in_imag[i] * ft_real[ft];
        }
    }

    // out[] includes both polarisations, at adjacent indices
    out[2*idx]   = -out_real / K;
    out[2*idx+1] = -out_imag / K;

    __syncthreads();
}


/**
 * Invert the PFB by applying a resynthesis filter, using GPU acceleration.
 *
 * @param[in]   data_buffer_fine The fine-channelised voltages
 * @param       file_no An index corresponding to the second of the data
 * @param       npointing The number of pointings
 * @param       nsamples The number of time samples per second of input data
 * @param       nchan The number of fine channels
 * @param       npol The number of instrumental polarisations
 * @param       g The ipfb filter struct
 * @param[out]  data_buffer_vdif The reconstructed coarse-channelised voltages
 *
 * This function expects `data_buffer_fine` to be a 1D array of complex
 * voltages with indices following the ordering:
 *
 * ```
 *   pointing, samp, chan, pol
 * ```
 *
 * Although `data_buffer_fine` potentially contains 2 seconds' worth of data,
 * this function only inverts 1 second. The appropriate second is worked out
 * using `file_no`: if it is even, the first half of `data_buffer_fine` is
 * used; if odd, the second half.
 *
 * The IPFB is also given ntaps samples from the end of the previous second,
 * so in total the IPFB is given nsamp+ntaps samples worth of data.
 *
 * The output of the inversion is packed back into `data_buffer_vdif`, a 1D
 * array whose ordering is as follows:
 *
 * ```
 *   time, pol, complexity
 * ```
 *
 * This ordering is suited for immediate output to the VDIF format.
 *
 * It is assumed that the inverse filter coefficients have already been loaded
 * to the GPU.
 */
void cu_invert_pfb( gpuDoubleComplex *data_buffer_fine, int file_no,
                        int npointing, int nsamples, int nchan, int npol,
                        struct gpu_ipfb_arrays *g, float *data_buffer_vdif )
{
    // Setup input values:
    // The starting sample index is "ntaps" places from the end of the second
    // half of detected_beam if the file number is even, and "ntaps" places
    // from the end of the first half of detected_beam if the file number is
    // odd.
    int start_s = (file_no % 2 == 0 ? 2*nsamples - g->ntaps : nsamples - g->ntaps);

    int p,s_in,s,ch,pol,i,j;
    for (p = 0; p < npointing; p++)
    for (s_in = 0; s_in < nsamples + g->ntaps; s_in++)
    {
        // s_in is the sample index in the IPFB buffers, in_real and in_imag
        // s    is the sample index in data_buffer_fine
        s = (start_s + s_in) % (2*nsamples);
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
            {
                // Calculate the index for the IPFB bufffers, in_real and in_imag;
                i = B_IDX(p,s_in,ch,pol,nsamples+g->ntaps,nchan,npol);

                // Calculate the index for data_buffer_fine
                j = B_IDX(p,s,ch,pol,nsamples,nchan,npol);
                
                // Copy the data into the IPFB buffers
                // For the first chunk of data, there is no overlap from the previous
                // second so the first ntaps samples will be set to zero
                if (file_no == 0 && s_in < g->ntaps)
                {
                    g->in_real[i] = 0.0;
                    g->in_imag[i] = 0.0;
                }
                else
                {
                    // Here, the swapping of real <--> imag was made to fix the long-standing
                    // "VDIF bug". It's unclear why we *have* to do this, but it seems to
                    // work. There is probably some other swapping that happens secretly elsewhere.
                    g->in_real[i] = gpuCimag( data_buffer_fine[j] );
                    g->in_imag[i] = gpuCreal( data_buffer_fine[j] );
                }
            }
        }
    }
    
    // Copy the data to the device
    (gpuMemcpy( g->d_in_real, g->in_real, g->in_size, gpuMemcpyHostToDevice ));
    (gpuMemcpy( g->d_in_imag, g->in_imag, g->in_size, gpuMemcpyHostToDevice ));
    
    // Call the kernel
    if (npointing > 1)
    {
        fprintf( stderr, "error: PFB inversion currently only supports a single pointing\n" );
        exit(EXIT_FAILURE);
    }
    ipfb_kernel<<<nsamples, nchan*npol>>>( g->d_in_real, g->d_in_imag,
                                             g->d_ft_real, g->d_ft_imag,
                                             g->ntaps, npol, g->d_out );
    ( gpuPeekAtLastError() );

    // Copy the result back into host memory
    (gpuMemcpy( data_buffer_vdif, g->d_out, g->out_size, gpuMemcpyDeviceToHost ));
}


/**
 * Load the inverse PFB filter coefficients to the GPU.
 *
 * This function loads the inverse filter coefficients and the twiddle factors
 * into GPU memory. If they were loaded separately (as floats), then the
 * multiplication of the filter coefficients and the twiddle factors will be
 * less precise than if a single array containing every combination of floats
 * and twiddle factors is calculated in doubles, and then demoted to floats.
 * Hence, this pre-calculation is done in this function before gpuMemcpy is
 * called.
 *
 * The result is 2x 1D arrays loaded onto the GPU (one for real, one for imag)
 * where the ith element is equal to
 * \f[
 *     \text{result}[i] = f[n] \times exp{2\pi jk/K},
 * \f]
 * where \f$n = i \% N\f$, \f$k = i / N\f$, \f$N\f$ is the filter size
 * (`fil_size`) and \f$K\f$ is the number of channels (`nchan`).
 */
void cu_load_ipfb_filter( pfb_filter *filter, struct gpu_ipfb_arrays *g )
{
    int ch, f, i;

    // Setup filter values:
    gpuDoubleComplex ft; // pre-calculated filter coeffs times twiddle factor
    gpuDoubleComplex cf; // temp variable for complex version of filter coeffs
    for (f = 0; f < filter->ncoeffs; f++)
    {
        cf = make_gpuDoubleComplex( filter->coeffs[f], 0.0 );
        for (ch = 0; ch < filter->nchans; ch++)
        {
            i = filter->ncoeffs*ch + f;
            ft = gpuCmul( filter->twiddles[ch], cf );
            g->ft_real[i] = gpuCreal( ft );
            g->ft_imag[i] = gpuCimag( ft );
        }
    }

    (gpuMemcpy( g->d_ft_real, g->ft_real, g->ft_size, gpuMemcpyHostToDevice ));
    (gpuMemcpy( g->d_ft_imag, g->ft_imag, g->ft_size, gpuMemcpyHostToDevice ));
}


/**
 * Allocate GPU memory needed for performing the inverse PFB.
 *
 * Free with free_ipfb().
 */
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
    (gpuMalloc( (void **)&g->d_in_real, g->in_size ));
    (gpuMalloc( (void **)&g->d_in_imag, g->in_size ));
    (gpuMalloc( (void **)&g->d_ft_real, g->ft_size ));
    (gpuMalloc( (void **)&g->d_ft_imag, g->ft_size ));

    (gpuMalloc( (void **)&g->d_out, g->out_size ));

    // Allocate memory for host copies of the same
    g->in_real = (float *)malloc( g->in_size );
    g->in_imag = (float *)malloc( g->in_size );
    g->ft_real = (float *)malloc( g->ft_size );
    g->ft_imag = (float *)malloc( g->ft_size );

}

/**
 * Free the GPU memory allocated in malloc_ipfb().
 */
void free_ipfb( struct gpu_ipfb_arrays *g )
{
    // Free memory on host and device
    free( g->in_real );
    free( g->in_imag );
    free( g->ft_real );
    free( g->ft_imag );
    gpuFree( g->d_in_real );
    gpuFree( g->d_in_imag );
    gpuFree( g->d_ft_real );
    gpuFree( g->d_ft_imag );
    gpuFree( g->d_out );
}
