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


__global__ void invj_the_data( uint8_t       *data,
                               cuDoubleComplex *J,
                               cuDoubleComplex *phi,
                               cuDoubleComplex *JDq,
                               cuDoubleComplex *JDp,
                               uint32_t      *polQ_idxs,
                               uint32_t      *polP_idxs,
                               int npol )
/* Layout for input arrays:
*   data [nsamples] [nchan] [ninputs]            -- see docs
*   J    [nants] [nchan] [npol] [npol]        -- jones matrix
*   incoh --true if outputing an incoherent beam
* Layout for output arrays:
*   JDq  [nsamples] [nchan] [nant]
*   JDp  [nsamples] [nchan] [nant]
*/
{
    // Translate GPU block/thread numbers into meaning->l names
    int c    = blockIdx.x;  /* The (c)hannel number */
    int nc   = gridDim.x;   /* The (n)umber of (c)hannels */
    int s    = blockIdx.y;  /* The (s)ample number */
    int nant = blockDim.x;  /* The (n)umber of (a)ntennas */

    int ant  = threadIdx.x; /* The (ant)enna number */

    int ni   = nant*npol;   /* The (n)umber of RF (i)nputs */

    int iQ   = polQ_idxs[ant]; /* The input index for the X pol for this antenna */
    int iP   = polP_idxs[ant]; /* The input index for the Y pol for this antenna */

    cuDoubleComplex Dq, Dp;
    // Convert input data to complex float
    Dq  = UCMPLX4_TO_CMPLX_FLT(data[v_IDX(s,c,iQ,nc,ni)]);
    Dp  = UCMPLX4_TO_CMPLX_FLT(data[v_IDX(s,c,iP,nc,ni)]);

    // Calculate the first step (J*D) of the coherent beam (B = J*phi*D)
    // JDq = Jqq*Dq + Jqp*Dp
    // JDp = Jpq*Dq + Jpy*Dp
    JDq[JD_IDX(s,c,ant,nc,nant)] = cuCadd( cuCmul( J[J_IDX(ant,c,0,0,nc,npol)], Dq ),
                                           cuCmul( J[J_IDX(ant,c,0,1,nc,npol)], Dp ) );
    JDp[JD_IDX(s,c,ant,nc,nant)] = cuCadd( cuCmul( J[J_IDX(ant,c,1,0,nc,npol)], Dq ),
                                           cuCmul( J[J_IDX(ant,c,1,1,nc,npol)], Dp ) );
#ifdef DEBUG
    if (c==0 && s==0 && ant==0)
    {
        printf( "Jinv = [jyq, jyp; jxq; jxp]\n"
                "     = [%lf%+lf*i, %lf%+lf*i; %lf%+lf*i, %lf%+lf*i]\n",
                cuCreal(J[J_IDX(ant,c,0,0,nc,npol)]), cuCimag(J[J_IDX(ant,c,0,0,nc,npol)]),
                cuCreal(J[J_IDX(ant,c,0,1,nc,npol)]), cuCimag(J[J_IDX(ant,c,0,1,nc,npol)]),
                cuCreal(J[J_IDX(ant,c,1,0,nc,npol)]), cuCimag(J[J_IDX(ant,c,1,0,nc,npol)]),
                cuCreal(J[J_IDX(ant,c,1,1,nc,npol)]), cuCimag(J[J_IDX(ant,c,1,1,nc,npol)]) );
        printf( "v    = [vq; vp] = [%.1lf%+.1lf*i; %.1lf%+.1lf*i]\n",
                cuCreal( Dq ), cuCimag( Dq ),
                cuCreal( Dp ), cuCimag( Dp ) );
    }
#endif
}

__global__ void beamform_kernel( cuDoubleComplex *JDq,
                                 cuDoubleComplex *JDp,
                                 cuDoubleComplex *phi,
                                 double invw,
                                 int p,
                                 int soffset,
                                 int nchunk,
                                 cuDoubleComplex *Bd,
                                 float *C,
                                 int npol )
/* Layout for input arrays:
*   JDq  [nsamples] [nchan] [nant]               -- calibrated voltages
*   JDp  [nsamples] [nchan] [nant]
*   phi  [nant    ] [nchan]                      -- weights array
*   invw                                         -- inverse atrix
* Layout of input options
*   p                                            -- pointing number
*   soffset                                      -- sample offset (10000/nchunk)
*   nchunk                                       -- number of chunks each second is split into
* Layout for output arrays:
*   Bd   [nsamples] [nchan]   [npol]             -- detected beam
*   C    [nsamples] [nstokes] [nchan]            -- coherent full stokes
*/
{
    // Translate GPU block/thread numbers into meaning->l names
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

    cuDoubleComplex *Bx  = (cuDoubleComplex *)(&arrays[1*nant]);
    cuDoubleComplex *By  = (cuDoubleComplex *)(&arrays[3*nant]);
    cuDoubleComplex *Nxx = (cuDoubleComplex *)(&arrays[5*nant]);
    cuDoubleComplex *Nxy = (cuDoubleComplex *)(&arrays[7*nant]);
    cuDoubleComplex *Nyy = (cuDoubleComplex *)(&arrays[9*nant]);

    // Calculate the beam and the noise floor

    /* Fix from Maceij regarding NaNs in output when running on Athena, 13 April 2018.
    Apparently the different compilers and architectures are treating what were
    unintialised variables very differently */

    Bx[ant]  = make_cuDoubleComplex( 0.0, 0.0 );
    By[ant]  = make_cuDoubleComplex( 0.0, 0.0 );

    Nxx[ant] = make_cuDoubleComplex( 0.0, 0.0 );
    Nxy[ant] = make_cuDoubleComplex( 0.0, 0.0 );
    //Nyx[ant] = make_cuDoubleComplex( 0.0, 0.0 );
    Nyy[ant] = make_cuDoubleComplex( 0.0, 0.0 );

    // Calculate beamform products for each antenna, and then add them together
    // Calculate the coherent beam (B = J*phi*D)
    Bx[ant] = cuCmul( phi[PHI_IDX(p,ant,c,nant,nc)], JDq[JD_IDX(s,c,ant,nc,nant)] );
    By[ant] = cuCmul( phi[PHI_IDX(p,ant,c,nant,nc)], JDp[JD_IDX(s,c,ant,nc,nant)] );

    Nxx[ant] = cuCmul( Bx[ant], cuConj(Bx[ant]) );
    Nxy[ant] = cuCmul( Bx[ant], cuConj(By[ant]) );
    //Nyx[ant] = cuCmul( By[ant], cuConj(Bx[ant]) ); // Not needed as it's degenerate with Nxy[]
    Nyy[ant] = cuCmul( By[ant], cuConj(By[ant]) );

    // Detect the coherent beam
    // A summation over an array is faster on a GPU if you add half on array
    // to its other half as than can be done in parallel. Then this is repeated
    // with half of the previous array until the array is down to 1.
    __syncthreads();
    for ( int h_ant = nant / 2; h_ant > 0; h_ant = h_ant / 2 )
    {
        if (ant < h_ant)
        {
            Bx[ant]  = cuCadd( Bx[ant],  Bx[ant  + h_ant] );
            By[ant]  = cuCadd( By[ant],  By[ant  + h_ant] );
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
        float bnXX = DETECT(Bx[0]) - cuCreal(Nxx[0]);
        float bnYY = DETECT(By[0]) - cuCreal(Nyy[0]);
        cuDoubleComplex bnXY = cuCsub( cuCmul( Bx[0], cuConj( By[0] ) ),
                                    Nxy[0] );

        // Stokes I, Q, U, V:
        C[C_IDX(p,s+soffset,0,c,ns,NSTOKES,nc)] = invw*(bnXX + bnYY);
        C[C_IDX(p,s+soffset,1,c,ns,NSTOKES,nc)] = invw*(bnXX - bnYY);
        C[C_IDX(p,s+soffset,2,c,ns,NSTOKES,nc)] =  2.0*invw*cuCreal( bnXY );
        C[C_IDX(p,s+soffset,3,c,ns,NSTOKES,nc)] = -2.0*invw*cuCimag( bnXY );

        // The beamformed products
        Bd[B_IDX(p,s+soffset,c,0,ns,nc,npol)] = Bx[0];
        Bd[B_IDX(p,s+soffset,c,1,ns,nc,npol)] = By[0];
    }
}

__global__ void renormalise_channels_kernel( float *I, int nstep, float *offsets, float *scales, uint8_t *Iscaled )
{
    // For just doing stokes I
    // One block
    // 128 threads; each thread will do one channel
    // (we have already summed over all ant)

    // For doing the C array (I,Q,U,V)
    // ... figure it out later.

    // Translate GPU block/thread numbers into meaningful names
    int chan    = threadIdx.x; /* The (c)hannel number */
    int nchan   = blockDim.x;  /* The total number of channels */

    int stokes  = threadIdx.y; /* The (stokes) parameter */
    int nstokes = blockDim.y;  /* Typically, is either 1 (just Stokes I) or 4 (Stokes IQUV) */

    int p       = blockIdx.x;  /* The (p)ointing number */

    float val, scale, offset;
    //float summed = 0.0;

    // Initialise min and max values to the first sample
    float min = I[C_IDX(p,0,stokes,chan,nstep,nstokes,nchan)];
    float max = I[C_IDX(p,0,stokes,chan,nstep,nstokes,nchan)];

    // Get the data statistics
    int i;
    for (i = 0; i < nstep; i++)
    {
        val = I[C_IDX(p,i,stokes,chan,nstep,nstokes,nchan)];
        min = (val < min ? val : min);
        max = (val > max ? val : max);
        //summed += fabsf(val);
    }

    // Rescale the incoherent beam to fit the available 8 bits
    scale  = (max - min) / 256.0; // (for 8-bit values)
    offset = min + 0.5*scale;

    for (i = 0; i < nstep; i++)
    {
        val = (I[C_IDX(p,i,stokes,chan,nstep,nstokes,nchan)] - offset) / scale;
        Iscaled[C_IDX(p,i,stokes,chan,nstep,nstokes,nchan)] = (uint8_t)(val + 0.5);
    }

    // Set the scales and offsets
    scales[p*nstokes*nchan + stokes*nchan + chan] = scale;
    offsets[p*nstokes*nchan + stokes*nchan + chan] = offset;

    // Originally, this function's purpose was to "flatten the bandpass";
    // however, it appears that all downstream pulsar software does this
    // routinely anyway. Along with the commented-out lines above relating
    // to the calculation of the mean, the following line was the original
    // "flattening" instruction, with "32" being a magic number to get
    // the range to fit nicely (NB but not perfectly!) into 8 bits.
    // As it is, Iscaled as calculated above, only affects the PSRFITS
    // output, while uncommenting the following line (along with the other
    // necessary lines above) will only affect the VDIF output.
    //
    // In particular, I currently believe that the PFB inversion should
    // in principle have higher fidelity if the bandpass is _not_ flattened
    // first, since the inverse filter is designed to be the complement of
    // the forward filter anyway. However, preliminary tests show that the
    // difference it marginal, even negligible. Uncomment at your own risk!

    //float mean = summed / nstep;
    //I[C_IDX(p,i,stokes,chan,nstep,nstokes,nchan)] *= 32.0/mean;
}



void cu_form_incoh_beam(
        uint8_t *data, uint8_t *d_data, size_t data_size, // The input data
        float *d_incoh, // The data, summed
        unsigned int nsample, int nchan, int ninput, // The data dimensions
        float *offsets, float *d_offsets, // data statistics: offsets
        float *scales, float *d_scales,   // data statistics: scales
        uint8_t *Iscaled, uint8_t *d_Iscaled, size_t Iscaled_size // The scaled answer
        )
/* Inputs:
*   data    = host array of 4bit+4bit complex numbers. For data order, refer to the
*             documentation.
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

void cu_form_beam( uint8_t *data, unsigned int sample_rate,
                   cuDoubleComplex *d_phi,
                   int file_no,
                   int npointing, int nants, int nchan,
                   int npol, double invw,
                   struct gpu_formbeam_arrays *g,
                   cuDoubleComplex ****detected_beam, float *coh,
                   cudaStream_t *streams, int nchunk,
                   mpi_psrfits *mpfs )
/* Inputs:
*   data    = array of 4bit+4bit complex numbers. For data order, refer to the
*             documentation.
*   sample_rate = The voltage sample rate, in Hz
*   d_phi   = geometric delays (device)  [pointing][nants][nchan]
*   J       = inverse Jones matrix.  [nants][nchan][npol][npol]
*   file_no = number of file we are processing, starting at 0.
*   nants     = number of antennas
*   nchan        = number of channels
*   npol         = number of polarisations (=2)
*   invw         = the reciprocal of the sum of the antenna weights
*   g            = struct containing pointers to various arrays on
*                  both host and device
*
* Outputs:
*   detected_beam = result of beamforming operation, summed over antennas
*                   [2*nsamples][nchan][npol]
*   coh           = result in Stokes parameters (minus noise floor)
*                   [nsamples][nstokes][nchan]
*/
{
    // Copy the data to the device
    gpuErrchk(cudaMemcpyAsync( g->d_J,    g->J, g->J_size,    cudaMemcpyHostToDevice ));

    // Divide the gpu calculation into multiple time chunks so there is enough room on the GPU
    int p;
    for (int ichunk = 0; ichunk < nchunk; ichunk++)
    {
        //int dataoffset = ichunk * g->data_size / sizeof(uint8_t);
        gpuErrchk(cudaMemcpyAsync( g->d_data,
                                   data + ichunk * g->data_size / sizeof(uint8_t),
                                   g->data_size, cudaMemcpyHostToDevice ));

        // Call the kernels
        // samples_chan(index=blockIdx.x  size=gridDim.x,
        //              index=blockIdx.y  size=gridDim.y)
        // stat_point  (index=threadIdx.x size=blockDim.x,
        //              index=threadIdx.y size=blockDim.y)
        //dim3 samples_chan(sample_rate, nchan);
        dim3 chan_samples( nchan, sample_rate / nchunk );
        dim3 stat( nants );

        // convert the data and multiply it by J
        invj_the_data<<<chan_samples, stat>>>( g->d_data, g->d_J, d_phi, g->d_JDq, g->d_JDp,
                                               g->d_polQ_idxs, g->d_polP_idxs,
                                               npol );

        // Send off a parallel CUDA stream for each pointing
        for (p = 0; p < npointing; p++ )
        {
            // Call the beamformer kernel
            // To see how the 11*STATION double arrays are used, go to tag 11NSTATION
            beamform_kernel<<<chan_samples, stat, 11*nants*sizeof(double), streams[p]>>>( g->d_JDq, g->d_JDp,
                            d_phi, invw,
                            p, ichunk*sample_rate/nchunk, nchunk,
                            g->d_Bd, g->d_coh, npol );

            gpuErrchk( cudaPeekAtLastError() );
        }
    }
    gpuErrchk( cudaDeviceSynchronize() );


    // Flatten the bandpass
    float *d_offsets, *offsets;
    float *d_scales, *scales;
    uint8_t *d_Cscaled, *Cscaled;

    size_t offsets_size = npointing*nchan*NSTOKES*sizeof(float);
    size_t scales_size  = npointing*nchan*NSTOKES*sizeof(float);
    size_t Cscaled_size = npointing*mpfs[0].coarse_chan_pf.sub.bytes_per_subint;

    gpuErrchk(cudaMalloc( (void **)&d_offsets, offsets_size ));
    gpuErrchk(cudaMalloc( (void **)&d_scales,  scales_size ));
    gpuErrchk(cudaMalloc( (void **)&d_Cscaled, Cscaled_size ));

    gpuErrchk(cudaMallocHost( (void **)&offsets, offsets_size ));
    gpuErrchk(cudaMallocHost( (void **)&scales,  scales_size ));
    gpuErrchk(cudaMallocHost( (void **)&Cscaled, Cscaled_size ));

    dim3 chan_stokes(nchan, NSTOKES);
    renormalise_channels_kernel<<<npointing, chan_stokes, 0, streams[0]>>>( g->d_coh, sample_rate, d_offsets, d_scales, d_Cscaled );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    gpuErrchk(cudaMemcpy( offsets, d_offsets, offsets_size, cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( scales,  d_scales,  scales_size,  cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( Cscaled, d_Cscaled, Cscaled_size, cudaMemcpyDeviceToHost ));

    for (p = 0; p < npointing; p++)
    {
        memcpy( mpfs[p].coarse_chan_pf.sub.dat_offsets, &(offsets[p*nchan*NSTOKES]), nchan*NSTOKES*sizeof(float) );
        memcpy( mpfs[p].coarse_chan_pf.sub.dat_scales, &(scales[p*nchan*NSTOKES]), nchan*NSTOKES*sizeof(float) );
        memcpy( mpfs[p].coarse_chan_pf.sub.data, &(Cscaled[p*sample_rate*nchan*NSTOKES]), sample_rate*nchan*NSTOKES );
    }

    gpuErrchk(cudaFreeHost( offsets ));
    gpuErrchk(cudaFreeHost( scales ));
    gpuErrchk(cudaFreeHost( Cscaled ));

    gpuErrchk(cudaFree( d_offsets ));
    gpuErrchk(cudaFree( d_scales ));
    gpuErrchk(cudaFree( d_Cscaled ));

    // Copy the results back into host memory
    gpuErrchk(cudaMemcpyAsync( g->Bd, g->d_Bd,    g->Bd_size,    cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpyAsync( coh,   g->d_coh,   g->coh_size,   cudaMemcpyDeviceToHost ));

    // Copy the data back from Bd back into the detected_beam array
    // Make sure we put it back into the correct half of the array, depending
    // on whether this is an even or odd second.
    int offset, i;
    offset = file_no % 2 * sample_rate;

    // TODO: turn detected_beam into a 1D array
    int ch, s, pol;
    for (p   = 0; p   < npointing  ; p++  )
    for (s   = 0; s   < sample_rate; s++  )
    for (ch  = 0; ch  < nchan      ; ch++ )
    for (pol = 0; pol < npol       ; pol++)
    {
        i = p  * (npol*nchan*sample_rate) +
            s  * (npol*nchan)                   +
            ch * (npol)                         +
            pol;

        detected_beam[p][s+offset][ch][pol] = g->Bd[i];
    }
}

void malloc_formbeam( struct gpu_formbeam_arrays *g, vcsbeam_metadata *vm,
                      int *nchunk, float gpu_mem_gb, int outpol_coh,
                      int npointing, logger *log )
{
    size_t data_base_size;
    size_t JD_base_size;

    int sample_rate = vm->sample_rate;
    int nants       = vm->obs_metadata->num_ants;
    int nchan       = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    int npol        = vm->obs_metadata->num_ant_pols; // (X,Y)

    // Calculate array sizes for host and device
    g->coh_size    = npointing * sample_rate * outpol_coh * nchan * sizeof(float);
    data_base_size = sample_rate * nants * nchan * npol * sizeof(uint8_t);
    //g->data_size  = sample_rate * nants * nchan * npol / nchunk * sizeof(uint8_t);
    g->Bd_size     = npointing * sample_rate * nchan * npol * sizeof(cuDoubleComplex);
    size_t phi_size = npointing * nants * nchan * sizeof(cuDoubleComplex);
    g->J_size      = nants * nchan * npol * npol * sizeof(cuDoubleComplex);
    JD_base_size   = sample_rate * nants * nchan * sizeof(cuDoubleComplex);
    //g->JD_size    = sample_rate * nants * nchan / nchunk * sizeof(cuDoubleComplex);
    
    size_t gpu_mem;
    if ( gpu_mem_gb == -1.0f ) {
        // Find total GPU memory
        struct cudaDeviceProp gpu_properties;
        cudaGetDeviceProperties( &gpu_properties, 0 );
        gpu_mem = gpu_properties.totalGlobalMem;
        gpu_mem_gb = (float)gpu_mem / (float)(1024*1024*1024);
    }
    else {
        gpu_mem = (size_t)(gpu_mem_gb * (float)(1024*1024*1024));
    }


    // Work out how many chunks to split a second into so there is enough memory on the gpu
    *nchunk = 0;
    size_t gpu_mem_used = pow(10, 15); // 1 PB
    while ( gpu_mem_used > gpu_mem ) 
    {
        *nchunk += 1;
        // Make sure the nchunk is divisable by the samples
        while ( sample_rate%*nchunk != 0 )
        {
            *nchunk += 1;
        }
        gpu_mem_used = (phi_size + g->J_size + g->Bd_size + data_base_size / *nchunk +
                        g->coh_size + 3*JD_base_size / *nchunk);
    }
    float gpu_mem_used_gb = (float)gpu_mem_used / (float)(1024*1024*1024);

    char log_message[128];

    sprintf( log_message, "Splitting each second into %d chunks", *nchunk );
    logger_timed_message( log, log_message );

    sprintf( log_message, "%6.3f GB out of the total %6.3f GPU memory allocated",
                     gpu_mem_used_gb, gpu_mem_gb );
    logger_timed_message( log, log_message );

    g->data_size = data_base_size / *nchunk;
    g->JD_size   = JD_base_size / *nchunk;


    // Allocate host memory
    //g->J  = (cuDoubleComplex *)malloc( g->J_size );
    //g->Bd = (cuDoubleComplex *)malloc( g->Bd_size );
    cudaMallocHost( &g->J, g->J_size );
    cudaCheckErrors("cudaMallocHost J fail");
    cudaMallocHost( &g->Bd, g->Bd_size );
    cudaCheckErrors("cudaMallocHost Bd fail");

    sprintf( log_message, "coh_size   %9.3f MB GPU mem", (float)g->coh_size  / (float)(1024*1024) );
    logger_timed_message( log, log_message );

    sprintf( log_message, "data_size  %9.3f MB GPU mem", (float)g->data_size / (float)(1024*1024) );
    logger_timed_message( log, log_message );

    sprintf( log_message, "Bd_size    %9.3f MB GPU mem", (float)g->Bd_size   / (float)(1024*1024) );
    logger_timed_message( log, log_message );

    sprintf( log_message, "phi_size   %9.3f MB GPU mem", (float)phi_size     / (float)(1024*1024) );
    logger_timed_message( log, log_message );

    sprintf( log_message, "J_size     %9.3f MB GPU mem", (float)g->J_size    / (float)(1024*1024) );
    logger_timed_message( log, log_message );

    sprintf( log_message, "JD_size    %9.3f MB GPU mem", (float)g->JD_size*3 / (float)(1024*1024) );
    logger_timed_message( log, log_message );


    // Allocate device memory
    gpuErrchk(cudaMalloc( (void **)&g->d_J,     g->J_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_JDq,   g->JD_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_JDp,   g->JD_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_Bd,    g->Bd_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_data,  g->data_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_coh,   g->coh_size ));

    // Allocate memory on both host and device for polX and polY idx arrays
    g->pol_idxs_size = nants * sizeof(uint32_t);
    cudaMallocHost( &g->polQ_idxs, g->pol_idxs_size );
    cudaMallocHost( &g->polP_idxs, g->pol_idxs_size );
    gpuErrchk(cudaMalloc( (void **)&g->d_polQ_idxs, g->pol_idxs_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_polP_idxs, g->pol_idxs_size ));
}

void cu_upload_pol_idx_lists( struct gpu_formbeam_arrays *g )
{
    gpuErrchk(cudaMemcpy( g->d_polQ_idxs, g->polQ_idxs, g->pol_idxs_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( g->d_polP_idxs, g->polP_idxs, g->pol_idxs_size, cudaMemcpyHostToDevice ));
}

void free_formbeam( struct gpu_formbeam_arrays *g )
{
    // Free memory on host and device
    cudaFreeHost( g->J );
    cudaFreeHost( g->Bd );
    cudaFreeHost( g->polQ_idxs );
    cudaFreeHost( g->polP_idxs );
    cudaFree( g->d_J );
    cudaFree( g->d_Bd );
    cudaFree( g->d_data );
    cudaFree( g->d_coh );
    cudaFree( g->d_polQ_idxs );
    cudaFree( g->d_polP_idxs );
}

float *create_pinned_data_buffer( size_t size )
{
    float *ptr;
    cudaMallocHost( &ptr, size );
    cudaCheckErrors("cudaMallocHost data_buffer fail");

    // Initialise to zeros
    memset( ptr, 0, size );

    return ptr;
}


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

void allocate_input_output_arrays( void **data, void **d_data, size_t size )
{
    cudaMallocHost( data, size );
    cudaCheckErrors( "cudaMallocHost() failed" );

    cudaMalloc( d_data, size );
    cudaCheckErrors( "cudaMalloc() failed" );
}

void free_input_output_arrays( void *data, void *d_data )
{
    cudaFreeHost( data );
    cudaCheckErrors( "cudaFreeHost() failed" );

    cudaFree( d_data );
    cudaCheckErrors( "cudaFree() failed" );
}


void flatten_bandpass(int nstep, int nchan, int npol, void *data)
{

    // nstep -> ? number of time steps (ie 10,000 per 1 second of data)
    // nchan -> number of (fine) channels (128)
    // npol -> number of polarisations (1 for icoh or 4 for coh)

    // magical mystery normalisation constant
    int new_var = 32;

    // purpose is to generate a mean value for each channel/polaridation

    int i=0, j=0;
    int p=0;

    float *data_ptr = (float *) data;
    float **band;


    band = (float **) calloc (npol, sizeof(float *));
    for (i=0;i<npol;i++) {
      band[i] = (float *) calloc(nchan, sizeof(float));
    }

    // initialise the band array
    for (p = 0;p<npol;p++) {
        for (j=0;j<nchan;j++){
            band[p][j] = 0.0;
        }
    }


    // accumulate abs(data) over all time samples and save into band
    data_ptr = (float *)data;
    for (i=0;i<nstep;i++) { // time steps
        for (p = 0;p<npol;p++) { // pols
            for (j=0;j<nchan;j++){ // channels
                band[p][j] += fabsf(*data_ptr);
                data_ptr++;
            }
        }

    }

    // calculate and apply the normalisation to the data
    data_ptr = (float *)data;
    for (i=0;i<nstep;i++) {
        for (p = 0;p<npol;p++) {
            for (j=0;j<nchan;j++){
                *data_ptr = (*data_ptr)/( (band[p][j]/nstep)/new_var );
                data_ptr++;
            }
        }

    }

    // free the memory
    for (i=0;i<npol;i++) {
        free(band[i]);
    }
    free(band);
}

