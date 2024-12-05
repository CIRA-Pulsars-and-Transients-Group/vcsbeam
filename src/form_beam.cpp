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

#include "gpu_fft.hpp"
#include "gpu_macros.h"


#include "vcsbeam.h"

/* Wrapper function for GPU/CUDA error handling.
 * 
 * @todo Remove this function and replace all calls to it with calls to CUDA
 *       error functions.
 */
/*inline void gpuAssert(gpuError_t code, const char *file, int line, bool abort=true)
{
    if (code != 0)
    {
        fprintf(stderr, "GPUAssert:: %s - %s (%d)\n", cudaGetErrorString(code), file, line);
        if (abort)
        {
            exit(code);
        }
    }
}*/

// define a macro for accessing gpuAssert
// #define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__, true);}


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
    gpuDoubleComplex v; // Complex voltage
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
 * @param p         The pointing number
 * @param soffset   An offset number of samples into `data`
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
                                 gpuDoubleComplex *J,
                                 gpuDoubleComplex *Jv_Q,
                                 gpuDoubleComplex *Jv_P,
                                 uint32_t      *polQ_idxs,
                                 uint32_t      *polP_idxs,
                                 int npol,
                                 int p,
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
    int ns   = gridDim.y;   /* The (n)umber of (s)amples (in a chunk) */
    int nant = blockDim.x;  /* The (n)umber of (a)ntennas */

    int ant  = threadIdx.x; /* The (ant)enna number */

    int ni   = nant*npol;   /* The (n)umber of RF (i)nputs */

    int iQ   = polQ_idxs[ant]; /* The input index for the Q pol for this antenna */
    int iP   = polP_idxs[ant]; /* The input index for the P pol for this antenna */

    gpuDoubleComplex vq, vp;
    // Convert input data to complex float
    if (datatype == VM_INT4)
    {
        uint8_t *v = (uint8_t *)data;
        vq = UCMPLX4_TO_CMPLX_FLT(v[v_IDX(s,c,iQ,nc,ni)]);
        vp = UCMPLX4_TO_CMPLX_FLT(v[v_IDX(s,c,iP,nc,ni)]);
    }
    else if (datatype == VM_DBL)
    {
        gpuDoubleComplex *v = (gpuDoubleComplex *)data;
        vq = v[v_IDX(s,c,iQ,nc,ni)];
        vp = v[v_IDX(s,c,iP,nc,ni)];
    }
    // else send an error message... yet to do

    // Calculate the first step (J*v) of the coherent beam
    // Jv_Q = Jqq*vq + Jqp*vp
    // Jv_P = Jpq*vq + Jpy*vp

    Jv_Q[Jv_IDX(p,s,c,ant,ns,nc,nant)] = gpuCadd( gpuCmul( J[J_IDX(p,ant,c,0,0,nant,nc,npol)], vq ),
                                                  gpuCmul( J[J_IDX(p,ant,c,0,1,nant,nc,npol)], vp ) );
    Jv_P[Jv_IDX(p,s,c,ant,ns,nc,nant)] = gpuCadd( gpuCmul( J[J_IDX(p,ant,c,1,0,nant,nc,npol)], vq ),
                                                  gpuCmul( J[J_IDX(p,ant,c,1,1,nant,nc,npol)], vp ) );

/*#ifdef DEBUG
    if (c==50 && s == 3 && ant==0)
    {
        for (int i = 0; i < 1; i++)
        {
            printf( "Jinv[%3d]      = [%5.3lf,%5.3lf, %5.3lf,%5.3lf]\n"
                    "                 [%5.3lf,%5.3lf, %5.3lf,%5.3lf]\n",
                    i,
                    gpuCreal(J[J_IDX(p,i,c,0,0,nant,nc,npol)]), gpuCimag(J[J_IDX(p,i,c,0,0,nant,nc,npol)]),
                    gpuCreal(J[J_IDX(p,i,c,0,1,nant,nc,npol)]), gpuCimag(J[J_IDX(p,i,c,0,1,nant,nc,npol)]),
                    gpuCreal(J[J_IDX(p,i,c,1,0,nant,nc,npol)]), gpuCimag(J[J_IDX(p,i,c,1,0,nant,nc,npol)]),
                    gpuCreal(J[J_IDX(p,i,c,1,1,nant,nc,npol)]), gpuCimag(J[J_IDX(p,i,c,1,1,nant,nc,npol)]) );
            gpuDoubleComplex *v = (gpuDoubleComplex *)data;
            printf( "v[Q=%3d,P=%3d] = [%5.3lf,%5.3lf]; [%5.3lf,%5.3lf]\n",
                    polQ_idxs[i], polP_idxs[i],
                    gpuCreal( v[v_IDX(s,c,polQ_idxs[i],nc,ni)] ), gpuCimag( v[v_IDX(s,c,polQ_idxs[i],nc,ni)] ),
                    gpuCreal( v[v_IDX(s,c,polP_idxs[i],nc,ni)] ), gpuCimag( v[v_IDX(s,c,polP_idxs[i],nc,ni)] )
                  );
            printf( "Jinv * v = [%5.3lf,%5.3lf]; [%5.3lf,%5.3lf]\n",
                    gpuCreal(Jv_Q[Jv_IDX(p,s,c,i,ns,nc,nant)]), gpuCimag(Jv_Q[Jv_IDX(p,s,c,i,ns,nc,nant)]), 
                    gpuCreal(Jv_P[Jv_IDX(p,s,c,i,ns,nc,nant)]), gpuCimag(Jv_P[Jv_IDX(p,s,c,i,ns,nc,nant)])
                  );
        }
    }
#endif */
}

/**
 * CUDA kernel for phasing up and summing the voltages over antenna.
 *
 * @param[in] nfine_chan Number of fine channels
 * @param[in] n_samples Number of time samples
 * @param[in] nant Number of antennas
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
 * @param nstokes  The number of stokes parameters to output
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
__global__ void vmBeamform_kernel(int nfine_chan,
                                 int n_samples, 
                                 int nant, 
                                 gpuDoubleComplex *Jv_Q,
                                 gpuDoubleComplex *Jv_P,
                                 gpuDoubleComplex *phi,
                                 double invw,
                                 int p,
                                 int soffset,
                                 int nchunk,
                                 gpuDoubleComplex *e,
                                 float *S,
                                 int npol,
                                 int nstokes )
{

    const unsigned int warp_id {threadIdx.x / warpSize};
    const unsigned int lane_id {threadIdx.x % warpSize};
    const unsigned int glb_warp_id {blockIdx.x * (blockDim.x / warpSize) + warp_id};
    __shared__ gpuDoubleComplex workspace[512 * 5];
    
    if(glb_warp_id >= n_samples * nfine_chan) return;
    

    // Translate GPU block/thread numbers into meaningful names
    int c    = glb_warp_id / n_samples;
    int nc   = nfine_chan;  /* The (n)umber of (c)hannels (=128) */
    int s    = glb_warp_id % n_samples; /* The (s)ample number */
    int ns   = n_samples;  /* The (n)umber of (s)amples (in a chunk)*/

    gpuDoubleComplex ex  = make_gpuDoubleComplex( 0.0, 0.0 );
    gpuDoubleComplex ey  = make_gpuDoubleComplex( 0.0, 0.0 );
    gpuDoubleComplex Nxx = make_gpuDoubleComplex( 0.0, 0.0 );
    gpuDoubleComplex Nxy = make_gpuDoubleComplex( 0.0, 0.0 );
    gpuDoubleComplex Nyy = make_gpuDoubleComplex( 0.0, 0.0 );
    // (Nyx is not needed as it's degenerate with Nxy)



    for(unsigned int ant {lane_id}; ant < nant; ant += warpSize){
        // Calculate beamform products for each antenna, and then add them together
        // Calculate the coherent beam (B = J*phi*D)
        gpuDoubleComplex ex_tmp = gpuCmul( phi[PHI_IDX(p,ant,c,nant,nc)], Jv_Q[Jv_IDX(p,s,c,ant,ns,nc,nant)] );
        gpuDoubleComplex ey_tmp = gpuCmul( phi[PHI_IDX(p,ant,c,nant,nc)], Jv_P[Jv_IDX(p,s,c,ant,ns,nc,nant)] );
        ex = gpuCadd(ex, ex_tmp);
        ey = gpuCadd(ey, ey_tmp);
        Nxx = gpuCadd(Nxx, gpuCmul( ex_tmp, gpuConj(ex_tmp)));
        Nxy = gpuCadd(Nxy, gpuCmul( ex_tmp, gpuConj(ey_tmp)));
        Nyy = gpuCadd(Nyy, gpuCmul( ey_tmp, gpuConj(ey_tmp)));
    }

    workspace[threadIdx.x * 5 + 0] = ex;
    workspace[threadIdx.x * 5 + 1] = ey;
    workspace[threadIdx.x * 5 + 2] = Nxx;
    workspace[threadIdx.x * 5 + 3] = Nxy;
    workspace[threadIdx.x * 5 + 4] = Nyy;

    
    for(unsigned int i = warpSize/2; i >= 1; i /= 2){
        if(lane_id < i){
            workspace[threadIdx.x * 5] = gpuCadd(workspace[threadIdx.x * 5], workspace[(threadIdx.x + i) * 5]);
            workspace[threadIdx.x * 5 + 1] = gpuCadd(workspace[threadIdx.x * 5 + 1], workspace[(threadIdx.x + i) * 5 + 1]);
            workspace[threadIdx.x * 5 + 2] = gpuCadd(workspace[threadIdx.x * 5 + 2], workspace[(threadIdx.x + i) * 5 + 2]);
            workspace[threadIdx.x * 5 + 3] = gpuCadd(workspace[threadIdx.x * 5 + 3], workspace[(threadIdx.x + i) * 5 + 3]);
            workspace[threadIdx.x * 5 + 4] = gpuCadd(workspace[threadIdx.x * 5 + 4], workspace[(threadIdx.x + i) * 5 + 4]);
        }
        #ifdef __NVCC__
        __syncwarp();
        #endif
    }
    // Form the stokes parameters for the coherent beam
    // Only doing it for ant 0 so that it only prints once
    if ( lane_id == 0 ) {
        float bnXX = DETECT(workspace[threadIdx.x * 5]) - gpuCreal(workspace[threadIdx.x * 5 + 2]);
        float bnYY = DETECT(workspace[threadIdx.x * 5 + 1]) - gpuCreal(workspace[threadIdx.x * 5 + 4]);
        gpuDoubleComplex bnXY = gpuCsub( gpuCmul( workspace[threadIdx.x * 5], gpuConj( workspace[threadIdx.x * 5 + 1] ) ),
                                    workspace[threadIdx.x * 5 + 3] );

        // Stokes I, Q, U, V:
        S[C_IDX(p,s+soffset,0,c,ns*nchunk,nstokes,nc)] = invw*(bnXX + bnYY);
        if ( nstokes == 4 )
        {
            S[C_IDX(p,s+soffset,1,c,ns*nchunk,nstokes,nc)] = invw*(bnXX - bnYY);
            S[C_IDX(p,s+soffset,2,c,ns*nchunk,nstokes,nc)] =  2.0*invw*gpuCreal( bnXY );
            S[C_IDX(p,s+soffset,3,c,ns*nchunk,nstokes,nc)] = -2.0*invw*gpuCimag( bnXY );
        }

        // The beamformed products
        e[B_IDX(p,s+soffset,c,0,ns*nchunk,nc,npol)] = workspace[threadIdx.x * 5 ];
        e[B_IDX(p,s+soffset,c,1,ns*nchunk,nc,npol)] = workspace[threadIdx.x * 5 + 1];
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
    scale  = (max - min) / 256.0; // (for 8-bit unsigned values)
    offset = min;

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
    (gpuMemcpyAsync( d_data, data, data_size, gpuMemcpyHostToDevice ));

    // Call the kernels
    dim3 chan_samples( nchan, nsample );

    // Call the incoherent beam kernel
    incoh_beam<<<chan_samples, ninput>>>( d_data, d_incoh );

    ( gpuPeekAtLastError() );
    ( gpuDeviceSynchronize() );

    dim3 chan_stokes(nchan, 1);
    int npointing = 1;
    renormalise_channels_kernel<<<npointing, chan_stokes>>>( d_incoh, nsample, d_offsets, d_scales, d_Iscaled );
    ( gpuPeekAtLastError() );

    // Copy the results back into host memory
    // (NB: Asynchronous copy here breaks the output)
    (gpuMemcpy( offsets, d_offsets, nchan*sizeof(float), gpuMemcpyDeviceToHost ));
    (gpuMemcpy( scales,  d_scales,  nchan*sizeof(float), gpuMemcpyDeviceToHost ));
    (gpuMemcpy( Iscaled, d_Iscaled, Iscaled_size,        gpuMemcpyDeviceToHost ));
}

/**
 * Computes \f${\bf J}^{-1} {\bf v}\f$.
 */
void vmApplyJChunk( vcsbeam_context *vm )
{
    dim3 chan_samples( vm->nfine_chan, vm->fine_sample_rate / vm->chunks_per_second );
    dim3 stat( vm->obs_metadata->num_ants );

    // J times v
    // Send off a parallel CUDA stream for each pointing
    int p;
    for (p = 0; p < vm->npointing; p++ )
    {
        vmApplyJ_kernel<<<chan_samples, stat, 0, vm->streams[p]>>>(
                vm->d_v,
                (gpuDoubleComplex *)vm->d_J,
                vm->d_Jv_Q,
                vm->d_Jv_P,
                vm->d_polQ_idxs,
                vm->d_polP_idxs,
                vm->obs_metadata->num_ant_pols,
                p,
                vm->datatype );
        gpuCheckLastError(); 
    }
    ( gpuDeviceSynchronize() );
}

/**
 * Performs the phasing up, averaging over antennas, and detection operations
 * on calibrated data.
 *
 * @todo Split the beamforming operations into separate steps/kernels.
 */
void vmBeamformChunk( vcsbeam_context *vm )
{
    /*
     * Cristian's implementation. Each warp (instead of an entire block) takes care of computing the
     beam (work item) for a frequency channel and time sample.
    */
    size_t total_work_items = vm->nfine_chan * (vm->fine_sample_rate / vm->chunks_per_second);
    const int warpSize = 64;
    const int nthreads = 512;
    const int warps_per_block = nthreads / warpSize;
    const int n_blocks = (total_work_items + warps_per_block - 1) / warps_per_block;

    // Get the "chunk" number
    int chunk = vm->chunk_to_load % vm->chunks_per_second;

    // Send off a parallel CUDA stream for each pointing
    int p;
    for (p = 0; p < vm->npointing; p++ )
    {
#ifdef DEBUG
        fprintf(stderr, "vm->npointing=%d  pointing=%d\n", vm->npointing, p);
        fprintf(stderr, "chan_samples=(%d,%d,%d) stat=(%d,%d,%d)\n", chan_samples.x, chan_samples.y, chan_samples.z, stat.x, stat.y, stat.z);
        fprintf(stderr, "I think the coarse channel numbers is: %d\n", vm->coarse_chan_idx);
#endif
        // Call the beamformer kernel
        vmBeamform_kernel<<<n_blocks, nthreads, 0, vm->streams[p]>>>(
                vm->nfine_chan, 
                (vm->fine_sample_rate / vm->chunks_per_second),
                vm->obs_metadata->num_ants,
                vm->d_Jv_Q,
                vm->d_Jv_P,
                vm->gdelays.d_phi,
                1.0/(double)vm->num_not_flagged,
                p,
                chunk*vm->fine_sample_rate/vm->chunks_per_second,
                vm->chunks_per_second,
                vm->d_e,
                (float *)vm->d_S,
                vm->obs_metadata->num_ant_pols,
                vm->out_nstokes );
        gpuCheckLastError();
    }
    ( gpuDeviceSynchronize() );
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
    ( gpuDeviceSynchronize() );

    // Unlock the buffer for reading
    // TODO: generalise this for arbitrary pipelines
    // SM: This custom "lock" implementation might prove futile in the face of asynchronous CUDA calls...
    vm->v->locked = false;
}

/**
 * Copies the beamformed voltages from GPU memory to CPU memory.
 */
void vmPullE( vcsbeam_context *vm )
{
    // Copy the results back into host memory
    gpuMemcpyAsync( vm->e, vm->d_e, vm->e_size_bytes, gpuMemcpyDeviceToHost );
}

/**
 * Copies the detected Stokes parameters from GPU memory to CPU memory.
 */
void vmPullS( vcsbeam_context *vm )
{
    gpuMemcpyAsync( vm->S, vm->d_S, vm->S_size_bytes, gpuMemcpyDeviceToHost );
}

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
    dim3 chan_stokes(vm->nfine_chan, vm->out_nstokes);
    renormalise_channels_kernel<<<vm->npointing, chan_stokes, 0, vm->streams[0]>>>( (float *)vm->d_S, vm->fine_sample_rate, vm->d_offsets, vm->d_scales, vm->d_Cscaled );
    ( gpuPeekAtLastError() );
    ( gpuDeviceSynchronize() );

    (gpuMemcpy( vm->offsets, vm->d_offsets, vm->offsets_size, gpuMemcpyDeviceToHost ));
    (gpuMemcpy( vm->scales,  vm->d_scales,  vm->scales_size,  gpuMemcpyDeviceToHost ));
    (gpuMemcpy( vm->Cscaled, vm->d_Cscaled, vm->Cscaled_size, gpuMemcpyDeviceToHost ));

    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        memcpy( mpfs[p].coarse_chan_pf.sub.dat_offsets, &(vm->offsets[p*vm->nfine_chan*vm->out_nstokes]), vm->nfine_chan*vm->out_nstokes*sizeof(float) );
        memcpy( mpfs[p].coarse_chan_pf.sub.dat_scales, &(vm->scales[p*vm->nfine_chan*vm->out_nstokes]), vm->nfine_chan*vm->out_nstokes*sizeof(float) );
        memcpy( mpfs[p].coarse_chan_pf.sub.data, &(vm->Cscaled[p*vm->fine_sample_rate*vm->nfine_chan*vm->out_nstokes]), vm->fine_sample_rate*vm->nfine_chan*vm->out_nstokes );
    }

}

/**
 * Copies the index arrays for antennas and polarisations from CPU memory to
 * GPU memory.
 */
void vmPushPolIdxLists( vcsbeam_context *vm )
{
    gpuMemcpy( vm->d_polQ_idxs, vm->polQ_idxs, vm->pol_idxs_size_bytes, gpuMemcpyHostToDevice );
    gpuMemcpy( vm->d_polP_idxs, vm->polP_idxs, vm->pol_idxs_size_bytes, gpuMemcpyHostToDevice );
}

/**
 * (Deprecated) Allocate memory on the GPU.
 *
 * @todo Remove the function create_pinned_data_buffer().
 */
float *create_pinned_data_buffer( size_t size )
{
    float *ptr;
    gpuMallocHost( &ptr, size );

    // Initialise to zeros
    memset( ptr, 0, size );

    return ptr;
}


/**
 * (Deprecated) Allocate memory on the CPU.
 *
 * @todo Remove the function create_detected_beam().
 */
gpuDoubleComplex ****create_detected_beam( int npointing, int nsamples, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, s, ch; // Loop variables
    gpuDoubleComplex ****array;

    array = (gpuDoubleComplex ****)malloc( npointing * sizeof(gpuDoubleComplex ***) );
    for (p = 0; p < npointing; p++)
    {
        array[p] = (gpuDoubleComplex ***)malloc( nsamples * sizeof(gpuDoubleComplex **) );

        for (s = 0; s < nsamples; s++)
        {
            array[p][s] = (gpuDoubleComplex **)malloc( nchan * sizeof(gpuDoubleComplex *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][s][ch] = (gpuDoubleComplex *)malloc( npol * sizeof(gpuDoubleComplex) );
        }
    }
    return array;
}


/*
* Create a 1-D host buffer containing the fine-channelised complex voltage samples which
* will be given to the IPFB.
*/
gpuDoubleComplex *create_data_buffer_fine( int npointing, int nsamples, int nchan, int npol )
{
    int buffer_size = npointing * nsamples * nchan * npol * sizeof(gpuDoubleComplex);
    
    // Allocate host memory for buffer
    gpuDoubleComplex *data_buffer_fine;
    data_buffer_fine = (gpuDoubleComplex *)malloc( buffer_size );

    // Initialise buffer to zeros
    memset( data_buffer_fine, 0, buffer_size );

    return data_buffer_fine;
}


/*
* Copy the fine-channelised beamformed voltages into a data buffer so that it can
* be given to the IPFB.
*/
void prepare_data_buffer_fine( gpuDoubleComplex *data_buffer_fine, vcsbeam_context *vm,
                    uintptr_t timestep_idx )
{
    // Get shortcut variables
    uintptr_t nchan = vm->nfine_chan;
    uintptr_t npol  = vm->obs_metadata->num_ant_pols; // = 2
    int file_no = timestep_idx % 2;

    // Copy the beamformed data from e into the data buffer
    // Make sure we put it back into the correct half of the array, depending
    // on whether this is an even or odd second.
    int offset = file_no % 2 * vm->fine_sample_rate;

    int p,s,ch,pol,i,j;
    for (p   = 0; p   < vm->npointing;        p++   )
    for (s   = 0; s   < vm->fine_sample_rate; s++   )
    for (ch  = 0; ch  < nchan;                ch++  )
    for (pol = 0; pol < npol;                 pol++ )
    {
        // Calculate index for e
        i = B_IDX(p,s,ch,pol,vm->fine_sample_rate,nchan,npol);

        // Calculate index for data_buffer_fine
        j = B_IDX(p,s+offset,ch,pol,vm->fine_sample_rate,nchan,npol);

        data_buffer_fine[j] = vm->e[i];
    }
}

/**
 * (Deprecated) Allocates memory on the CPU and GPU simultaneously.
 *
 * @todo Remove the function allocate_input_output_arrays().
 */
void allocate_input_output_arrays( void **data, void **d_data, size_t size )
{
    gpuMallocHost( data, size );

    gpuMalloc( d_data, size );
}

/**
 * (Deprecated) Frees memory on the CPU and GPU simultaneously.
 *
 * @todo Remove the function free_input_output_arrays().
 */
void free_input_output_arrays( void *data, void *d_data )
{
    gpuHostFree( data );

    gpuFree( d_data );
}
