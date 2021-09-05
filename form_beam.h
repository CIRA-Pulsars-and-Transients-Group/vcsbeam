/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <cuda_runtime.h>
#include <cuComplex.h>

#ifndef FORM_BEAM_H
#define FORM_BEAM_H

#include "jones.h"
#include "performance.h"

/* structure for managing data arrays to be allocated on both host and device */
struct gpu_formbeam_arrays
{
    size_t coh_size;
    size_t incoh_size;
    size_t data_size;
    size_t Bd_size;
    size_t J_size;
    size_t JD_size;
    size_t pol_idxs_size;
    cuDoubleComplex *J, *d_J;
    cuDoubleComplex *Bd, *d_Bd;
    cuDoubleComplex *JDx, *d_JDx;
    cuDoubleComplex *JDy, *d_JDy;
    uint8_t *d_data;
    float   *d_coh;
    float   *d_incoh;
    float   *d_Ia;
    uint32_t *polX_idxs, *d_polX_idxs;
    uint32_t *polY_idxs, *d_polY_idxs;
};


void malloc_formbeam( struct gpu_formbeam_arrays *g, unsigned int sample_rate,
                      int nstation, int nchan, int npol, int *nchunk, float gpu_mem_gb, int outpol_coh,
                      int outpol_incoh, int npointing, logger *log );
void free_formbeam( struct gpu_formbeam_arrays *g );



/* Converting from 4+4 complex to full-blown complex doubles */

#define REAL_NIBBLE_TO_UINT8(X)  ((X) & 0xf)
#define IMAG_NIBBLE_TO_UINT8(X)  (((X) >> 4) & 0xf)
#define UINT8_TO_INT(X)          ((X) >= 0x8 ? (signed int)(X) - 0x10 : (signed int)(X))
#define RE_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))))
#define IM_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X))))
#define UCMPLX4_TO_CMPLX_FLT(X)  (make_cuDoubleComplex((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))), \
                                         (float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X)))))
#define DETECT(X)                (cuCreal(cuCmul(X,cuConj(X))))



void cu_form_incoh_beam(
        uint8_t *data, uint8_t *d_data, size_t data_size,
        float *d_incoh,
        unsigned int nsample, int nchan, int ninput,
        float *offsets, float *d_offsets,
        float *scales, float *d_scales,
        uint8_t *Iscaled, uint8_t *d_Iscaled, size_t Iscaled_size
        );


void cu_form_beam( uint8_t *data, unsigned int sample_rate, cuDoubleComplex *d_phi,
                   cuDoubleComplex ****J, int file_no, 
                   int npointing, int nstation, int nchan,
                   int npol, double invw, struct gpu_formbeam_arrays *g,
                   cuDoubleComplex ****detected_beam, float *coh,
                   cudaStream_t *streams, int nchunk  );

float *create_pinned_data_buffer_psrfits( size_t size );
        
float *create_pinned_data_buffer_vdif( size_t size );

void cu_upload_pol_idx_lists( struct gpu_formbeam_arrays *g );

cuDoubleComplex ****create_invJi( int nstation, int nchan, int pol );
void              destroy_invJi( cuDoubleComplex ****array, int nstation, int nchan, int npol );

cuDoubleComplex ****create_detected_beam( int npointing, int nsamples, int nchan, int npol );
void              destroy_detected_beam( cuDoubleComplex ****array, int npointing,
                                         int nsamples, int nchan );

void allocate_input_output_arrays( void **data, void **d_data, size_t size );
void free_input_output_arrays( void *data, void *d_data );

void flatten_bandpass( int nstep, int nchan, int npol, void *data);


#endif
