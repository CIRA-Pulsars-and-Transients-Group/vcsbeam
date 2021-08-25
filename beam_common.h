/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef BEAM_COMMON_H
#define BEAM_COMMON_H

#include <inttypes.h>
#include "mwa_hyperbeam.h"
#include "calibration.h"
#include <cuComplex.h>
#include <mwalib.h>

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



struct beam_geom {
    double mean_ra;
    double mean_dec;
    double az;
    double el;
    double lmst;
    double fracmjd;
    double intmjd;
    double unit_N;
    double unit_E;
    double unit_H;
};


/* Running get_jones from within make_beam */
void get_jones(
        // an array of pointings [pointing][ra/dec][characters]
        int                    npointing, // number of pointings
        MetafitsMetadata      *osb_metadata,
        int                    coarse_chan_idx,
        struct                 calibration *cal,
        cuDoubleComplex     ***D,
        cuDoubleComplex       *B,
        cuDoubleComplex    ****invJi                   // output
);

void create_antenna_lists( MetafitsMetadata *metafits_metadata, uint32_t *polX_idxs, uint32_t *polY_idxs );

void int8_to_uint8(int n, int shift, char * to_convert);
void float2int8_trunc(float *f, int n, float min, float max, int8_t *i);
void float_to_unit8(float * in, int n, int8_t *out);

void dec2hms( char *out, double in, int sflag );
void utc2mjd( char *, double *, double * ); // "2000-01-01T00:00:00" --> MJD_int + MJD_fraction
void mjd2lst( double, double * );

double parse_dec( char* ); // "01:23:45.67" --> Dec in degrees
double parse_ra( char* );  // "01:23:45.67" --> RA  in degrees

/**** MATRIX OPERATIONS ****/

void cp2x2(cuDoubleComplex *Min, cuDoubleComplex *Mout);
void inv2x2(cuDoubleComplex *Min, cuDoubleComplex *Mout);
void inv2x2d(double *Min, double *Mout);
void inv2x2S(cuDoubleComplex *Min, cuDoubleComplex **Mout);
void mult2x2d(cuDoubleComplex *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout);
void mult2x2d_RxC(double *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout);
void conj2x2(cuDoubleComplex *M, cuDoubleComplex *Mout);
double norm2x2(cuDoubleComplex *M, cuDoubleComplex *Mout);

void parse_pointing_file( const char *filename, double **ras_hours, double **decs_degs, unsigned int *npointings );

#endif
