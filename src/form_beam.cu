/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include <cuComplex.h>
#include <cuda_runtime.h>

#include "vcsbeam.h"

/* Wrapper function for GPU/CUDA error handling.
 * 
 * @todo Remove this function and replace all calls to it with calls to CUDA
 *       error functions.
 */
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
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


/**
 * CUDA kernel for computing an incoherent beam.
 *
 * @param[in] data The voltage data, \f$v\f$, with layout \f$N_t \times N_f \times N_i\f$.
 * @param[out] incoh The detected (Stokes I) powers, \f$I\f$, with layout \f$N_t \times N_f\f$.
 *
 * The incoherent beam is the expression
 * \f[
 * I_{t,f} = \sum_i v_{t,f,i}^\dagger v_{t,f,i}.
 * \f]
 *
 * The expected thread configuration is
 * \f$\langle\langle\langle(N_f, N_t), N_i\rangle\rangle\rangle.\f$
 *
 * Note that if the voltages were arranged into Jones vectors, the above could
 * also be expressed in the more familiar form
 * \f[
 * I_{t,f} = \sum_a {\bf v}_{t,f,a}^\dagger {\bf v}_{t,f,a}.
 * \f]
 */
__global__ void incoh_beam( uint8_t *data, float *incoh )
/* <<< (nchan,nsample), ninput >>>
 */
{
    // Translate GPU block/thread numbers into meaningful names
    int c    = blockIdx.x;  /* The (c)hannel number */
    int nc   = gridDim.x;   /* The (n)umber of (c)hannels */
    int s    = blockIdx.y;  /* The (s)ample number */
    int ni   = blockDim.x;  /* The (n)umber of RF (i)nputs */

    int i  = threadIdx.x; /* The (ant)enna number */

    int idx = I_IDX(s, c, nc); /* Index into incoh */

    if (i == 0)
        incoh[idx] = 0.0;
    __syncthreads();

    // Convert input data to complex double
    cuDoubleComplex v; // Complex voltage
    uint8_t sample = data[v_IDX(s,c,i,nc,ni)];
    v = UCMPLX4_TO_CMPLX_FLT(sample);

    // Detect the sample ("detect" = calculate power = magnitude squared)
    // and add it to the others from this thread
    atomicAdd( &incoh[idx], DETECT(v) );
    __syncthreads();
}


/**
 * CUDA kernel for multiplying Jones matrices to Jones vectors.
 *
 * @param[in] data The voltage data, \f$v\f$,
 *                 with layout \f$N_t \times N_f \times N_i\f$
 * @param[in] J The Jones matrices, \f${\bf J}^{-1}\f$,
 *              with layout \f$N_a \times N_f \times N_p \times N_p\f$
 * @param[out] Jv_Q The Q polarisation of the product \f${\bf J}^{-1}{\bf v}\f$,
 *             with layout \f$N_t \times N_f \times N_a\f$
 * @param[out] Jv_P The P polarisation of the product \f${\bf J}^{-1}{\bf v}\f$,
 *             with layout \f$N_t \times N_f \times N_a\f$
 * @param polQ_idxs And array of the indices \f$i\f$ for the Q polarisations of
 *                  the antennas
 * @param polP_idxs And array of the indices \f$i\f$ for the P polarisations of
 *                  the antennas
 * @param npol      \f$N_p\f$
 * @param datatype Either `VM_INT4` (if `data` contain 4+4-bit complex integers)
 *                 or `VM_DBL` (if `data` contain complex doubles).
 *
 * Although this kernel is quite general, in the sense that it could be used
 * to multiply any Jones matrices to any Jones vectors, it is used in particular
 * for multiplying the Jones matrices \f${\bf J}^{-1}\f$ to the voltage data
 * \f${\bf v}\f$:
 * \f[
 * \tilde{\bf e}_{t,f,a} = {\bf J}^{-1}_{a,f}{\bf v}_{t,f,a}.
 * \f]
 *
 * The expected thread configuration is
 * \f$\langle\langle\langle(N_f, N_t), N_a\rangle\rangle\rangle.\f$
 */
__global__ void vmApplyJ_kernel( void            *data,
                                 cuDoubleComplex *J,
                                 cuDoubleComplex *Jv_Q,
                                 cuDoubleComplex *Jv_P,
                                 uint32_t      *polQ_idxs,
                                 uint32_t      *polP_idxs,
                                 int npol,
                                 vcsbeam_datatype datatype )
/* Layout for input arrays:
*   data [nsamples] [nchan] [ninputs]            -- see docs
*   J    [nants] [nchan] [npol] [npol]        -- jones matrix
*   incoh --true if outputing an incoherent beam
* Layout for output arrays:
*   Jv_Q  [nsamples] [nchan] [nant]
*   Jv_P  [nsamples] [nchan] [nant]
*/
{
    // Translate GPU block/thread numbers into meaningful names
    int c    = blockIdx.x;  /* The (c)hannel number */
    int nc   = gridDim.x;   /* The (n)umber of (c)hannels */
    int s    = blockIdx.y;  /* The (s)ample number */
    int nant = blockDim.x;  /* The (n)umber of (a)ntennas */

    int ant  = threadIdx.x; /* The (ant)enna number */

    int ni   = nant*npol;   /* The (n)umber of RF (i)nputs */

    int iQ   = polQ_idxs[ant]; /* The input index for the Q pol for this antenna */
    int iP   = polP_idxs[ant]; /* The input index for the P pol for this antenna */

    cuDoubleComplex vq, vp;
    // Convert input data to complex float
    if (datatype == VM_INT4)
    {
        uint8_t *v = (uint8_t *)data;
        vq = UCMPLX4_TO_CMPLX_FLT(v[v_IDX(s,c,iQ,nc,ni)]);
        vp = UCMPLX4_TO_CMPLX_FLT(v[v_IDX(s,c,iP,nc,ni)]);
    }
    else if (datatype == VM_DBL)
    {
        cuDoubleComplex *v = (cuDoubleComplex *)data;
        vq = v[v_IDX(s,c,iQ,nc,ni)];
        vp = v[v_IDX(s,c,iP,nc,ni)];
    }
    // else send an error message... yet to do

    // Calculate the first step (J*v) of the coherent beam
    // Jv_Q = Jqq*vq + Jqp*vp
    // Jv_P = Jpq*vq + Jpy*vp
    Jv_Q[Jv_IDX(s,c,ant,nc,nant)] = cuCadd( cuCmul( J[J_IDX(ant,c,0,0,nc,npol)], vq ),
                                           cuCmul( J[J_IDX(ant,c,0,1,nc,npol)], vp ) );
    Jv_P[Jv_IDX(s,c,ant,nc,nant)] = cuCadd( cuCmul( J[J_IDX(ant,c,1,0,nc,npol)], vq ),
                                           cuCmul( J[J_IDX(ant,c,1,1,nc,npol)], vp ) );
#ifdef DEBUG
    if (c==50 && s==0 && ant==0)
    {
        printf( "Jinv = [jyq, jyp; jxq; jxp]\n"
                "     = [%lf%+lf*i, %lf%+lf*i; %lf%+lf*i, %lf%+lf*i]\n",
                cuCreal(J[J_IDX(ant,c,0,0,nc,npol)]), cuCimag(J[J_IDX(ant,c,0,0,nc,npol)]),
                cuCreal(J[J_IDX(ant,c,0,1,nc,npol)]), cuCimag(J[J_IDX(ant,c,0,1,nc,npol)]),
                cuCreal(J[J_IDX(ant,c,1,0,nc,npol)]), cuCimag(J[J_IDX(ant,c,1,0,nc,npol)]),
                cuCreal(J[J_IDX(ant,c,1,1,nc,npol)]), cuCimag(J[J_IDX(ant,c,1,1,nc,npol)]) );
        printf( "v    = [vq; vp] = [%.1lf%+.1lf*i; %.1lf%+.1lf*i]\n",
                cuCreal( vq ), cuCimag( vq ),
                cuCreal( vp ), cuCimag( vp ) );
    }
#endif
}

/**
 * CUDA kernel for phasing up and summing the voltages over antenna.
 *
 * @param[in] Jv_Q The Q polarisation of the product \f${\bf J}^{-1}{\bf v}\f$,
 *                 with layout \f$N_t \times N_f \times N_a\f$
 * @param[in] Jv_P The P polarisation of the product \f${\bf J}^{-1}{\bf v}\f$,
 *                 with layout \f$N_t \times N_f \times N_a\f$
 * @param[in] phi  The delay phase, \f$\varphi\f$,
 *                 with layout \f$N_a \times N_f\f$
 * @param invw     The reciprocal of the number of non-flagged antennas
 * @param p        The pointing index
 * @param soffset  An offset number of samples into `e` for where to put the
 *                 the answer
 * @param nchunk   The number of chunks (divisions of a second's worth of data)
 * @param[out] e   The recovered electric field, \f${\bf e}\f$,
 *                 with layout \f$N_t \times N_f \times N_p\f$
 * @param[out] S   The recovered Stokes parameters,
 *                 with layout \f$N_t \times N_s \times N_f\f$
 * @param npol     \f$N_p\f$
 *
 * This kernel performs the phasing up and the summing over antennas part of
 * the beamforming operation (see [Beamforming](@ref beamforming)):
 * \f[
 *     {\bf e}_{t,f} = \frac{1}{N_a} \sum_a e^{i\varphi} \tilde{\bf e}_{t,f,a}.
 * \f]
 *
 * It also computes the Stokes parameters, \f$S = [I, Q, U, V]\f$ (with the
 * autocorrelations removed).
 *
 * The expected thread configuration is
 * \f$\langle\langle\langle(N_f, N_t), N_a\rangle\rangle\rangle.\f$
 */
__global__ void vmBeamform_kernel( cuDoubleComplex *Jv_Q,
                                 cuDoubleComplex *Jv_P,
                                 cuDoubleComplex *phi,
                                 double invw,
                                 int p,
                                 int soffset,
                                 int nchunk,
                                 cuDoubleComplex *e,
                                 float *S,
                                 int npol )
{
    // Translate GPU block/thread numbers into meaningful names
    int c    = blockIdx.x;  /* The (c)hannel number */
    int nc   = gridDim.x;   /* The (n)umber of (c)hannels (=128) */
    int s    = blockIdx.y;  /* The (s)ample number */
    int ns   = gridDim.y*nchunk;   /* The (n)umber of (s)amples (=10000)*/

    int ant  = threadIdx.x; /* The (ant)enna number */
    int nant = blockDim.x;  /* The (n)_umber of (ant)ennas */

    /*// GPU profiling
    clock_t start, stop;
    double setup_t, detect_t, sum_t, stokes_t;
    if ((p == 0) && (ant == 0) && (c == 0) && (s == 0)) start = clock();*/

    // Organise dynamically allocated shared arrays (see tag 11NSTATION for kernel call)
    extern __shared__ double arrays[];

    cuDoubleComplex *ex  = (cuDoubleComplex *)(&arrays[1*nant]);
    cuDoubleComplex *ey  = (cuDoubleComplex *)(&arrays[3*nant]);
    cuDoubleComplex *Nxx = (cuDoubleComplex *)(&arrays[5*nant]);
    cuDoubleComplex *Nxy = (cuDoubleComplex *)(&arrays[7*nant]);
    cuDoubleComplex *Nyy = (cuDoubleComplex *)(&arrays[9*nant]);

    // Calculate the beam and the noise floor

    /* Fix from Maceij regarding NaNs in output when running on Athena, 13 April 2018.
    Apparently the different compilers and architectures are treating what were
    unintialised variables very differently */

    ex[ant]  = make_cuDoubleComplex( 0.0, 0.0 );
    ey[ant]  = make_cuDoubleComplex( 0.0, 0.0 );

    Nxx[ant] = make_cuDoubleComplex( 0.0, 0.0 );
    Nxy[ant] = make_cuDoubleComplex( 0.0, 0.0 );
    //Nyx[ant] = make_cuDoubleComplex( 0.0, 0.0 );
    Nyy[ant] = make_cuDoubleComplex( 0.0, 0.0 );

    // Calculate beamform products for each antenna, and then add them together
    // Calculate the coherent beam (B = J*phi*D)
    ex[ant] = cuCmul( phi[PHI_IDX(p,ant,c,nant,nc)], Jv_Q[Jv_IDX(s,c,ant,nc,nant)] );
    ey[ant] = cuCmul( phi[PHI_IDX(p,ant,c,nant,nc)], Jv_P[Jv_IDX(s,c,ant,nc,nant)] );

    Nxx[ant] = cuCmul( ex[ant], cuConj(ex[ant]) );
    Nxy[ant] = cuCmul( ex[ant], cuConj(ey[ant]) );
    //Nyx[ant] = cuCmul( ey[ant], cuConj(ex[ant]) ); // Not needed as it's degenerate with Nxy[]
    Nyy[ant] = cuCmul( ey[ant], cuConj(ey[ant]) );

#ifdef DEBUG
    if (c==50 && s==0 && ant==0)
    {
        printf( "e    = [ex; ey] = [%lf%+lf*i; %lf%+lf*i]\n",
                cuCreal( ex[ant] ), cuCimag( ex[ant] ),
                cuCreal( ey[ant] ), cuCimag( ey[ant] ) );
        printf( "phi  = %lf%+lf*i\n",
                cuCreal(phi[PHI_IDX(p,ant,c,nant,nc)]),
                cuCimag(phi[PHI_IDX(p,ant,c,nant,nc)]) );
    }
#endif

    // Detect the coherent beam
    // A summation over an array is faster on a GPU if you add half on array
    // to its other half as than can be done in parallel. Then this is repeated
    // with half of the previous array until the array is down to 1.
    __syncthreads();
    for ( int h_ant = nant / 2; h_ant > 0; h_ant = h_ant / 2 )
    {
        if (ant < h_ant)
        {
            ex[ant]  = cuCadd( ex[ant],  ex[ant  + h_ant] );
            ey[ant]  = cuCadd( ey[ant],  ey[ant  + h_ant] );
            Nxx[ant] = cuCadd( Nxx[ant], Nxx[ant + h_ant] );
            Nxy[ant] = cuCadd( Nxy[ant], Nxy[ant + h_ant] );
            //Nyx[ant]=cuCadd( Nyx[ant], Nyx[ant + h_ant] );
            Nyy[ant] = cuCadd( Nyy[ant], Nyy[ant + h_ant] );
        }
        // below makes no difference so removed
        // else return;
        __syncthreads();
    }

    // Form the stokes parameters for the coherent beam
    // Only doing it for ant 0 so that it only prints once
    if ( ant == 0 )
    {
        float bnXX = DETECT(ex[0]) - cuCreal(Nxx[0]);
        float bnYY = DETECT(ey[0]) - cuCreal(Nyy[0]);
        cuDoubleComplex bnXY = cuCsub( cuCmul( ex[0], cuConj( ey[0] ) ),
                                    Nxy[0] );

        // Stokes I, Q, U, V:
        S[C_IDX(p,s+soffset,0,c,ns,NSTOKES,nc)] = invw*(bnXX + bnYY);
        S[C_IDX(p,s+soffset,1,c,ns,NSTOKES,nc)] = invw*(bnXX - bnYY);
        S[C_IDX(p,s+soffset,2,c,ns,NSTOKES,nc)] =  2.0*invw*cuCreal( bnXY );
        S[C_IDX(p,s+soffset,3,c,ns,NSTOKES,nc)] = -2.0*invw*cuCimag( bnXY );

        // The beamformed products
        e[B_IDX(p,s+soffset,c,0,ns,nc,npol)] = ex[0];
        e[B_IDX(p,s+soffset,c,1,ns,nc,npol)] = ey[0];
    }
}

/**
 * CUDA kernel for normalising Stokes parameters
 *
 * @param[in]  S       The original Stokes parameters,
 *                     with layout \f$N_t \times N_s \times N_f\f$
 * @param      nstep   \f$N_t\f$
 * @param[out] offsets The amount of offset needed to recover the original
 *                     values from the normalised ones
 * @param[out] scales  The scaling needed to recover the original values
 *                     from the normalised ones
 * @param[out] Sscaled The normalised Stokes parameters
 *
 * This kernel shifts and normalises the Stokes parameters so that they fit
 * into 8-bits integers without clipping (e.g. for output into the PSRFITS
 * format). Each frequency is normalised independently, with the scales and
 * offsets needed to recover the original values for that channel being
 * recorded as well.
 *
 * If \f${\bf S}\f$ is the array of values to be normalised, then the
 * normalisation is
 * \f[
 *     \hat{\bf S} = \frac{{\bf S} - \text{offset}}{\text{scale}},
 * \f]
 * where
 * \f{align*}{
 *     \text{scale} &= \frac{S_\text{max} - S_\text{min}}{256} \\
 *     \text{offset} &= S_\text{min} + 0.5 \times \text{scale}.
 * \f}
 *
 * The expected thread configuration is
 * \f$\langle\langle\langle N_b,(N_f, N_s)\rangle\rangle\rangle.\f$
 *
 * @todo Optimise the renormalisation kernel (e.g. by removing the for loop
 *       over timesteps.
 */
__global__ void renormalise_channels_kernel( float *S, int nstep, float *offsets, float *scales, uint8_t *Sscaled )
{
    // Translate GPU block/thread numbers into meaningful names
    int chan    = threadIdx.x; /* The (c)hannel number */
    int nchan   = blockDim.x;  /* The total number of channels */

    int stokes  = threadIdx.y; /* The (stokes) parameter */
    int nstokes = blockDim.y;  /* Typically, is either 1 (just Stokes I) or 4 (Stokes IQUV) */

    int p       = blockIdx.x;  /* The (p)ointing number */

    float val, scale, offset;
    //float summed = 0.0;

    // Initialise min and max values to the first sample
    float min = S[C_IDX(p,0,stokes,chan,nstep,nstokes,nchan)];
    float max = S[C_IDX(p,0,stokes,chan,nstep,nstokes,nchan)];

    // Get the data statistics
    int i;
    for (i = 0; i < nstep; i++)
    {
        val = S[C_IDX(p,i,stokes,chan,nstep,nstokes,nchan)];
        min = (val < min ? val : min);
        max = (val > max ? val : max);
        //summed += fabsf(val);
    }

    // Rescale the incoherent beam to fit the available 8 bits
    scale  = (max - min) / 256.0; // (for 8-bit values)
    offset = min + 0.5*scale;

    for (i = 0; i < nstep; i++)
    {
        val = (S[C_IDX(p,i,stokes,chan,nstep,nstokes,nchan)] - offset) / scale;
        Sscaled[C_IDX(p,i,stokes,chan,nstep,nstokes,nchan)] = (uint8_t)(val + 0.5);
    }

    // Set the scales and offsets
    scales[p*nstokes*nchan + stokes*nchan + chan] = scale;
    offsets[p*nstokes*nchan + stokes*nchan + chan] = offset;
}


/**
 * Form an incoherent beam.
 *
 * Forms an incoherent beam, detects it, and prepares it for writing to
 * PSRFITS.
 *
 * @todo Either generalise the vmBeamformChunk() so that it can also produce
 *       incoherent beams (and therefore do away with cu_form_incoh_beam(), or
 *       keep cu_form_incoh_beam() and convert it into a bona fide "vm" style
 *       fuction.
 */
void cu_form_incoh_beam(
        uint8_t *data, uint8_t *d_data, size_t data_size, // The input data
        float *d_incoh, // The data, summed
        unsigned int nsample, int nchan, int ninput, // The data dimensions
        float *offsets, float *d_offsets, // data statistics: offsets
        float *scales, float *d_scales,   // data statistics: scales
        uint8_t *Iscaled, uint8_t *d_Iscaled, size_t Iscaled_size // The scaled answer
        )
/* Inputs:
*   data    = host array of 4bit+4bit complex numbers.
*   d_data  = device array of the above
*   data_size    = size in bytes of data and d_data
*   nsample = number of samples
*   ninput  = number of RF inputs
*   nchan   = number of channels
*
* Outputs:
*   incoh  = result (Stokes I)   [nsamples][nchan]
*/
{
    // Copy the data to the device
    gpuErrchk(cudaMemcpyAsync( d_data, data, data_size, cudaMemcpyHostToDevice ));

    // Call the kernels
    dim3 chan_samples( nchan, nsample );

    // Call the incoherent beam kernel
    incoh_beam<<<chan_samples, ninput>>>( d_data, d_incoh );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    dim3 chan_stokes(nchan, 1);
    int npointing = 1;
    renormalise_channels_kernel<<<npointing, chan_stokes>>>( d_incoh, nsample, d_offsets, d_scales, d_Iscaled );
    gpuErrchk( cudaPeekAtLastError() );

    // Copy the results back into host memory
    // (NB: Asynchronous copy here breaks the output)
    gpuErrchk(cudaMemcpy( offsets, d_offsets, nchan*sizeof(float), cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( scales,  d_scales,  nchan*sizeof(float), cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( Iscaled, d_Iscaled, Iscaled_size,        cudaMemcpyDeviceToHost ));
}

/**
 * Computes \f${\bf J}^{-1} {\bf v}\f$.
 */
void vmApplyJChunk( vcsbeam_context *vm )
{
    dim3 chan_samples( vm->nfine_chan, vm->fine_sample_rate / vm->chunks_per_second );
    dim3 stat( vm->obs_metadata->num_ants );

    // J times v
    vmApplyJ_kernel<<<chan_samples, stat>>>(
            vm->d_v,
            (cuDoubleComplex *)vm->d_J,
            vm->d_Jv_Q,
            vm->d_Jv_P,
            vm->d_polQ_idxs,
            vm->d_polP_idxs,
            vm->obs_metadata->num_ant_pols,
            vm->datatype );
    cudaCheckErrors( "vmApplyJChunk: vmApplyJ_kernel failed" );
}

/**
 * Performs the phasing up, averaging over antennas, and detection operations
 * on calibrated data.
 *
 * @todo Split the beamforming operations into separate steps/kernels.
 */
void vmBeamformChunk( vcsbeam_context *vm )
{
    uintptr_t shared_array_size = 11 * vm->obs_metadata->num_ants * sizeof(double);
    // (To see how the 11*STATION double arrays are used, go to this code tag: 11NSTATION)

    dim3 chan_samples( vm->nfine_chan, vm->fine_sample_rate / vm->chunks_per_second );
    dim3 stat( vm->obs_metadata->num_ants );

    // Get the "chunk" number
    int chunk = vm->chunk_to_load % vm->chunks_per_second;
    // Send off a parallel CUDA stream for each pointing
    int p;
    for (p = 0; p < vm->npointing; p++ )
    {
        // Call the beamformer kernel
        vmBeamform_kernel<<<chan_samples, stat, shared_array_size, vm->streams[p]>>>(
                vm->d_Jv_Q,
                vm->d_Jv_P,
                vm->gdelays.d_phi,
                1.0/(double)vm->num_not_flagged,
                p,
                chunk*vm->fine_sample_rate/vm->chunks_per_second,
                vm->chunks_per_second,
                vm->d_e,
                (float *)vm->d_S,
                vm->obs_metadata->num_ant_pols );

        cudaCheckErrors( "vmBeamformChunk: vmBeamform_kernel failed" );
    }
    gpuErrchk( cudaDeviceSynchronize() );

}

/**
 * Performs all beamforming steps for 1 second's worth of data.
 */
void vmBeamformSecond( vcsbeam_context *vm )
{
    // Processing a second's worth of "chunks"
    int chunk;
    for (chunk = 0; chunk < vm->chunks_per_second; chunk++)
    {
        if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
        {
            vmPushChunk( vm );
        }
        else
        {
            // Do the forward PFB (see pfb.cu for better commented copy)
            vmUploadForwardPFBChunk( vm );

            logger_start_stopwatch( vm->log, "pfb", chunk == 0 ); // (report only on first round)

            vmWOLAChunk( vm );
            vmFPGARoundingChunk( vm );
            vmFFTChunk( vm );
            vmPackChunk( vm );

            logger_stop_stopwatch( vm->log, "pfb" );
        }

        logger_start_stopwatch( vm->log, "calc", chunk == 0 ); // (report only on first round)

        vmApplyJChunk( vm );
        vmBeamformChunk( vm );

        logger_stop_stopwatch( vm->log, "calc" );

        vm->chunk_to_load++;
    }
    gpuErrchk( cudaDeviceSynchronize() );

    // Unlock the buffer for reading
    // TODO: generalise this for arbitrary pipelines
    vm->v->locked = false;
}

/**
 * Copies the beamformed voltages from GPU memory to CPU memory.
 */
void vmPullE( vcsbeam_context *vm )
{
    // Copy the results back into host memory
    cudaMemcpyAsync( vm->e, vm->d_e, vm->e_size_bytes, cudaMemcpyDeviceToHost );
    cudaCheckErrors( "vmPullE: cudaMemcpyAsync failed" );
}

/**
 * Copies the detected Stokes parameters from GPU memory to CPU memory.
 */
void vmPullS( vcsbeam_context *vm )
{
    cudaMemcpyAsync( vm->S, vm->d_S, vm->S_size_bytes, cudaMemcpyDeviceToHost );
    cudaCheckErrors( "vmPullE: cudaMemcpyAsync failed" );
}

/**
 * Copies the beamformed voltages into a different array to prepare for
 * inverting the PFB.
 *
 * @todo Remove prepare_detected_beam() and use a "host buffer" instead.
 */
void prepare_detected_beam( cuDoubleComplex ****detected_beam,
        vcsbeam_context *vm )
{
    // Get shortcut variables
    uintptr_t nchan  = vm->nfine_chan;
    uintptr_t npol   = vm->obs_metadata->num_ant_pols; // = 2
    int file_no = vm->chunk_to_load / vm->chunks_per_second;

    // Copy the data back from e back into the detected_beam array
    // Make sure we put it back into the correct half of the array, depending
    // on whether this is an even or odd second.
    int offset, i;
    offset = file_no % 2 * vm->fine_sample_rate;

    // TODO: turn detected_beam into a 1D array
    int p, ch, s, pol;
    for (p   = 0; p   < vm->npointing  ; p++  )
    for (s   = 0; s   < vm->fine_sample_rate; s++  )
    for (ch  = 0; ch  < nchan      ; ch++ )
    for (pol = 0; pol < npol       ; pol++)
    {
        i = p  * (npol*nchan*vm->fine_sample_rate) +
            s  * (npol*nchan)                   +
            ch * (npol)                         +
            pol;

        detected_beam[p][s+offset][ch][pol] = vm->e[i];
    }
}

#ifdef HAVE_PSRFITS

/**
 * Renormalises the detected Stokes parameters and copies them into PSRFITS
 * structs, ready for frequency splicing.
 *
 * @param vm The VCSBeam context struct
 * @param mpfs The MPI PSRFITS struct that manages the splicing operation.
 */
void vmSendSToFits( vcsbeam_context *vm, mpi_psrfits *mpfs )
{
    // Flatten the bandpass
    dim3 chan_stokes(vm->nfine_chan, NSTOKES);
    renormalise_channels_kernel<<<vm->npointing, chan_stokes, 0, vm->streams[0]>>>( (float *)vm->d_S, vm->fine_sample_rate, vm->d_offsets, vm->d_scales, vm->d_Cscaled );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    gpuErrchk(cudaMemcpy( vm->offsets, vm->d_offsets, vm->offsets_size, cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( vm->scales,  vm->d_scales,  vm->scales_size,  cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( vm->Cscaled, vm->d_Cscaled, vm->Cscaled_size, cudaMemcpyDeviceToHost ));

    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        memcpy( mpfs[p].coarse_chan_pf.sub.dat_offsets, &(vm->offsets[p*vm->nfine_chan*NSTOKES]), vm->nfine_chan*NSTOKES*sizeof(float) );
        memcpy( mpfs[p].coarse_chan_pf.sub.dat_scales, &(vm->scales[p*vm->nfine_chan*NSTOKES]), vm->nfine_chan*NSTOKES*sizeof(float) );
        memcpy( mpfs[p].coarse_chan_pf.sub.data, &(vm->Cscaled[p*vm->fine_sample_rate*vm->nfine_chan*NSTOKES]), vm->fine_sample_rate*vm->nfine_chan*NSTOKES );
    }

}

#endif

/**
 * Copies the index arrays for antennas and polarisations from CPU memory to
 * GPU memory.
 */
void vmPushPolIdxLists( vcsbeam_context *vm )
{
    cudaMemcpy( vm->d_polQ_idxs, vm->polQ_idxs, vm->pol_idxs_size_bytes, cudaMemcpyHostToDevice );
    cudaCheckErrors( "vmMemcpyPolIdxLists: cudaMemcpy(polQ_idxs) failed" );
    cudaMemcpy( vm->d_polP_idxs, vm->polP_idxs, vm->pol_idxs_size_bytes, cudaMemcpyHostToDevice );
    cudaCheckErrors( "vmMemcpyPolIdxLists: cudaMemcpy(polP_idxs) failed" );
}

/**
 * (Deprecated) Allocate memory on the GPU.
 *
 * @todo Remove the function create_pinned_data_buffer().
 */
float *create_pinned_data_buffer( size_t size )
{
    float *ptr;
    cudaMallocHost( &ptr, size );
    cudaCheckErrors("cudaMallocHost data_buffer fail");

    // Initialise to zeros
    memset( ptr, 0, size );

    return ptr;
}


/**
 * (Deprecated) Allocate memory on the CPU.
 *
 * @todo Remove the function create_detected_beam().
 */
cuDoubleComplex ****create_detected_beam( int npointing, int nsamples, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, s, ch; // Loop variables
    cuDoubleComplex ****array;

    array = (cuDoubleComplex ****)malloc( npointing * sizeof(cuDoubleComplex ***) );
    for (p = 0; p < npointing; p++)
    {
        array[p] = (cuDoubleComplex ***)malloc( nsamples * sizeof(cuDoubleComplex **) );

        for (s = 0; s < nsamples; s++)
        {
            array[p][s] = (cuDoubleComplex **)malloc( nchan * sizeof(cuDoubleComplex *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][s][ch] = (cuDoubleComplex *)malloc( npol * sizeof(cuDoubleComplex) );
        }
    }
    return array;
}


/**
 * (Deprecated) Frees memory on the CPU.
 *
 * @todo Remove the function destroy_detected_beam().
 */
void destroy_detected_beam( cuDoubleComplex ****array, int npointing, int nsamples, int nchan )
{
    int p, s, ch;
    for (p = 0; p < npointing; p++)
    {
        for (s = 0; s < nsamples; s++)
        {
            for (ch = 0; ch < nchan; ch++)
                free( array[p][s][ch] );

            free( array[p][s] );
        }

        free( array[p] );
    }

    free( array );
}

/**
 * (Deprecated) Allocates memory on the CPU and GPU simultaneously.
 *
 * @todo Remove the function allocate_input_output_arrays().
 */
void allocate_input_output_arrays( void **data, void **d_data, size_t size )
{
    cudaMallocHost( data, size );
    cudaCheckErrors( "cudaMallocHost() failed" );

    cudaMalloc( d_data, size );
    cudaCheckErrors( "cudaMalloc() failed" );
}

/**
 * (Deprecated) Frees memory on the CPU and GPU simultaneously.
 *
 * @todo Remove the function free_input_output_arrays().
 */
void free_input_output_arrays( void *data, void *d_data )
{
    cudaFreeHost( data );
    cudaCheckErrors( "cudaFreeHost() failed" );

    cudaFree( d_data );
    cudaCheckErrors( "cudaFree() failed" );
}
