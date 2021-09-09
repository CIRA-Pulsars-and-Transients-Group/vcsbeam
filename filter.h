/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef FILTER_H
#define FILTER_H

#include <cuComplex.h>

#define REAL_COEFFS  0
#define CPLX_COEFFS  1

typedef enum filter_type_t
{
    ANALYSIS_FILTER,
    SYNTHESIS_FILTER
} filter_type;

typedef struct pfb_filter_t
{
    double          *coeffs;
    int              ncoeffs;
    int              ntaps;
    int              nchans; // = size/ntaps
    cuDoubleComplex *twiddles; // twiddle factors
    filter_type      type;
} pfb_filter;

#include <cuComplex.h>

cuDoubleComplex *roots_of_unity( int N );

pfb_filter *load_filter_coefficients( char *filtername, filter_type type, int nchans );
void free_pfb_filter( pfb_filter *filter );

#endif
