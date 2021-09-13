/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __JONES_H__
#define __JONES_H__

#include <inttypes.h>
#include "mwa_hyperbeam.h"
#include "calibration.h"
#include <cuComplex.h>
#include <mwalib.h>

#define NSTOKES  4

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

/* Calculating array indices for various flavours of Jones matrix arrays */

#define v_IDX(s,c,i,nc,ni)    ((s) * (ni)*(nc) + \
                               (c) * (ni)      + \
                               (i))

#define VSPVB  64000  /* "Vcsmwax Samples Per Voltage Block"
                         (i.e. the size of the TIME dimension in a voltage block */
#define vMWAX_IDX(s,i,ni) ((s)*(ni) + (i)*VSPVB + (s%VSPVB))

#define J_IDX(a,c,p1,p2,nc,npol)   ((a)  * (npol*npol*(nc))      + \
                                    (c)  * (npol*npol)           + \
                                    (p1) * (npol)                + \
                                    (p2))

#define JD_IDX(s,c,a,nc,na)  ((s) * (na)*(nc) + \
                              (c) * (na)      + \
                              (a))

#define B_IDX(p,s,c,pol,ns,nc,npol) ((p)  * (npol)*(nc)*(ns)   + \
                                     (s)  * (npol)*(nc)        + \
                                     (c)  * (npol)             + \
                                     (pol))

#define C_IDX(p,s,st,c,ns,nst,nc)  ((p)  * ((nc)*(nst)*(ns)) + \
                                    (s)  * ((nc)*(nst))      + \
                                    (st) *  (nc)               + \
                                    (c))

#define I_IDX(s,c,nc)          ((s)*(nc) + (c))



void create_antenna_lists( MetafitsMetadata *metafits_metadata, uint32_t *polX_idxs, uint32_t *polY_idxs );

void get_jones(
        // an array of pointings [pointing][ra/dec][characters]
        int                    npointing, // number of pointings
        MetafitsMetadata      *osb_metadata,
        int                    coarse_chan_idx,
        struct                 calibration *cal,
        cuDoubleComplex       *D,
        cuDoubleComplex       *B,
        cuDoubleComplex       *invJi // <-- output
);

/**** MATRIX OPERATIONS ****/

void cp2x2( cuDoubleComplex *Min, cuDoubleComplex *Mout );
void inv2x2( cuDoubleComplex *Min, cuDoubleComplex *Mout );
void inv2x2d( double *Min, double *Mout );
void inv2x2S( cuDoubleComplex *Min, cuDoubleComplex *Mout );
void mult2x2d( cuDoubleComplex *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout );
void mult2x2d_RxC( double *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout );
void conj2x2( cuDoubleComplex *M, cuDoubleComplex *Mout );
double norm2x2( cuDoubleComplex *M, cuDoubleComplex *Mout );
void calc_hermitian( cuDoubleComplex *M, cuDoubleComplex *H );
void calc_coherency_matrix( cuDoubleComplex *M, cuDoubleComplex *c );

#endif
