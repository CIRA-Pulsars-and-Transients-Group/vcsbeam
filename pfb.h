/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __PFB_H__
#define __PFB_H__

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>

#include <mwalib.h>

#include "filter.h"

/* The final step in the forward PFB version that emulates the FPGA
 * implementation packs the data into (4+4)-bit complex samples. The
 * following macros collectively achieve this.
 */
#define CLIP(x,max)        ((x) < -(max)   ? -(max)   : \
                            (x) >  (max)-1 ?  (max)-1 : (x))
#define INT_TO_UINT8(x)    ((x) < 0 ? CLIP(x,8) + 0x10 : CLIP(x,8))
#define DEMOTE(x)  (INT_TO_UINT8((int)round(x)))
#define PACK_NIBBLES(r,i)  ((DEMOTE(i) << 4) + DEMOTE(r))

/*******************************
 * FORWARD (ANALYSIS) FINE PFB *
 *******************************/

typedef struct forward_pfb_t
{
    char2            *htr_data;               // The input data, as obtained via mwalib from coarse channelised data
    char2            *htr_data_extended;      // For any extra "spillover" data beyond htr_data
    char2            *d_htr_data;             // All of the data, from both htr_data and htr_data_extended, on the device

    uint8_t          *vcs_data;               // The output data, fine channelised and packed into the VCS recombined format
    uint8_t          *d_vcs_data;             // Same as above, on device

    size_t            htr_size;               // The size (in bytes) of htr_data
    size_t            htr_extended_size;      // The size (in bytes) of htr_data_extended
    size_t            vcs_size;               // The size (in bytes) of vcs_data

    cuFloatComplex   *d_weighted_overlap_add; // A "temporary" array on the device for mid-calculation product

    int              *filter_coeffs;          // The filter to be applied
    int              *d_filter_coeffs;        // As above, on the device

    int               nspectra;               // The number of spectra to generate
    int               M;                      // The "stride" of the PFB (setting M=K means critically sampled)
    int               K;                      // The number of channels
    int               I;                      // The number of RF inputs
    int               P;                      // The number of taps

    cufftHandle       plan;                   // The cuFFT plan for performing the FFT part of the forward PFB
} forward_pfb;

//forward_pfb *init_forward_pfb( MetafitsMetadata *obs_metadata );
void cu_forward_pfb_fpga_version( forward_pfb *fpfb, bool copy_result_to_host );

/**********************************
 * BACKWARDS (SYNTHESIS) FINE PFB *
 **********************************/

struct gpu_ipfb_arrays
{
    int ntaps;
    int in_size;
    int ft_size;
    int out_size;
    float *in_real,   *in_imag;
    float *ft_real,   *ft_imag;
    float *d_in_real, *d_in_imag;
    float *d_ft_real, *d_ft_imag;
    float *d_out;
};

void cu_invert_pfb( cuDoubleComplex ****detected_beam, int file_no,
                        int npointing, int nsamples, int nchan, int npol, int sizeof_buffer,
                        struct gpu_ipfb_arrays *g, float *data_buffer_uvdif );

void cu_load_ipfb_filter( pfb_filter *filter, struct gpu_ipfb_arrays *g );

void malloc_ipfb( struct gpu_ipfb_arrays *g, pfb_filter *filter, int nsamples,
        int npol, int npointing );

void free_ipfb( struct gpu_ipfb_arrays *g );

#endif
