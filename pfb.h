/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __PFB_H__
#define __PFB_H__

#include <cuda_runtime.h>
#include <cuComplex.h>

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
