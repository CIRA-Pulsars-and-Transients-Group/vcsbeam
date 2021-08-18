/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuComplex.h>
#include "beam_common.h"
#include "psrfits.h"
#include "mwa_hyperbeam.h"

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
    data_ptr = data;
    for (i=0;i<nstep;i++) { // time steps
        for (p = 0;p<npol;p++) { // pols
            for (j=0;j<nchan;j++){ // channels
                band[p][j] += fabsf(*data_ptr);
                data_ptr++;
            }
        }

    }

    // calculate and apply the normalisation to the data
    data_ptr = data;
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

void int8_to_uint8(int n, int shift, char * to_convert) {
    int j;
    int scratch;
    int8_t with_sign;

    for (j = 0; j < n; j++) {
        with_sign = (int8_t) *to_convert;
        scratch = with_sign + shift;
        *to_convert = (uint8_t) scratch;
        to_convert++;
    }
}


void float2int8_trunc(float *f, int n, float min, float max, int8_t *i)
{
    int j;
    for (j = 0; j < n; j++) {
        f[j] = (f[j] > max) ? (max) : f[j];
        f[j] = (f[j] < min) ? (min) : f[j];
        i[j] = (int8_t) rint(f[j]);

    }
}

void float_to_unit8(float * in, int n, int8_t *out)
{
    int j;
    float min = -128.0; // -126.0 and -128.0 give the same result on test data
    float max = 127.0;
    // use a temp var so we don't modify the input data
    float scratch;
    for (j = 0; j < n; j++) {
        // TODO: count the number of samples that were clipped, store that and put it in the psrfits header
        // the if branching and ternary updates seem to be equivalent execution time
        if (in[j]> max) {
            scratch = max;
        } else if (in[j] < min) {
            scratch = min;
        } else {
            scratch = in[j];
        }
//        scratch = (in[j] > max) ? (max) : in[j];
//        scratch = (in[j] < min) ? (min) : scratch;
        out[j] = (uint8_t)( (int8_t)rint(scratch) + 128);
    }

}

void cp2x2(cuDoubleComplex *Min, cuDoubleComplex *Mout)
{
    Mout[0] = Min[0];
    Mout[1] = Min[1];
    Mout[2] = Min[2];
    Mout[3] = Min[3];
}


cuDoubleComplex reciprocal_complex( cuDoubleComplex z )
{
    double scale = 1.0/(z.x*z.x + z.y*z.y);
    return make_cuDoubleComplex( scale*z.x, -scale*z.y );
}

cuDoubleComplex negate_complex( cuDoubleComplex z )
{
    return make_cuDoubleComplex( -z.x, -z.y );
}

void inv2x2(cuDoubleComplex *Min, cuDoubleComplex *Mout)
{
    cuDoubleComplex m00 = Min[0];
    cuDoubleComplex m01 = Min[1];
    cuDoubleComplex m10 = Min[2];
    cuDoubleComplex m11 = Min[3];

    cuDoubleComplex m1 = cuCmul( m00, m11 );
    cuDoubleComplex m2 = cuCmul( m01, m10 );

    cuDoubleComplex det = cuCsub( m1, m2 );
    cuDoubleComplex inv_det = reciprocal_complex( det );

    Mout[0] = cuCmul(       inv_det,  m11 );
    Mout[1] = cuCmul( negate_complex(inv_det), m01 );
    Mout[2] = cuCmul( negate_complex(inv_det), m10 );
    Mout[3] = cuCmul(       inv_det,  m00 );
}

void inv2x2d(double *Min, double *Mout)
{
    double m00 = Min[0];
    double m01 = Min[1];
    double m10 = Min[2];
    double m11 = Min[3];

    double m1 = m00 * m11;
    double m2 = m01 * m10;

    double det = m1 - m2;
    double inv_det = 1/det;

    Mout[0] =  inv_det * m11;
    Mout[1] = -inv_det * m01;
    Mout[2] = -inv_det * m10;
    Mout[3] =  inv_det * m00;
}


void inv2x2S(cuDoubleComplex *Min, cuDoubleComplex **Mout)
// Same as inv2x2(), but the output is a 2x2 2D array, instead of a 4-element
// 1D array
{
    cuDoubleComplex m1 = cuCmul( Min[0], Min[3] );
    cuDoubleComplex m2 = cuCmul( Min[1], Min[2] );
    cuDoubleComplex det = cuCsub( m1, m2 );
    cuDoubleComplex inv_det = reciprocal_complex( det );
    Mout[0][0] = cuCmul(       inv_det,  Min[3] );
    Mout[0][1] = cuCmul( negate_complex(inv_det), Min[1] );
    Mout[1][0] = cuCmul( negate_complex(inv_det), Min[2] );
    Mout[1][1] = cuCmul(       inv_det,  Min[0] );
}


void mult2x2d(cuDoubleComplex *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout)
{
    cuDoubleComplex m00 = cuCmul( M1[0], M2[0] );
    cuDoubleComplex m12 = cuCmul( M1[1], M2[2] );
    cuDoubleComplex m01 = cuCmul( M1[0], M2[1] );
    cuDoubleComplex m13 = cuCmul( M1[1], M2[3] );
    cuDoubleComplex m20 = cuCmul( M1[2], M2[0] );
    cuDoubleComplex m32 = cuCmul( M1[3], M2[2] );
    cuDoubleComplex m21 = cuCmul( M1[2], M2[1] );
    cuDoubleComplex m33 = cuCmul( M1[3], M2[3] );
    Mout[0] = cuCadd( m00, m12 );
    Mout[1] = cuCadd( m01, m13 );
    Mout[2] = cuCadd( m20, m32 );
    Mout[3] = cuCadd( m21, m33 );
}

void mult2x2d_RxC(double *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout)
/* Mout = M1 x M2
 */
{
    cuDoubleComplex m00 = make_cuDoubleComplex( M1[0]*cuCreal(M2[0]), M1[0]*cuCimag(M2[0]) );
    cuDoubleComplex m12 = make_cuDoubleComplex( M1[1]*cuCreal(M2[2]), M1[1]*cuCimag(M2[2]) );
    cuDoubleComplex m01 = make_cuDoubleComplex( M1[0]*cuCreal(M2[1]), M1[0]*cuCimag(M2[1]) );
    cuDoubleComplex m13 = make_cuDoubleComplex( M1[1]*cuCreal(M2[3]), M1[1]*cuCimag(M2[3]) );
    cuDoubleComplex m20 = make_cuDoubleComplex( M1[2]*cuCreal(M2[0]), M1[2]*cuCimag(M2[0]) );
    cuDoubleComplex m32 = make_cuDoubleComplex( M1[3]*cuCreal(M2[2]), M1[3]*cuCimag(M2[2]) );
    cuDoubleComplex m21 = make_cuDoubleComplex( M1[2]*cuCreal(M2[1]), M1[2]*cuCimag(M2[1]) );
    cuDoubleComplex m33 = make_cuDoubleComplex( M1[3]*cuCreal(M2[3]), M1[3]*cuCimag(M2[3]) );
    Mout[0] = cuCadd( m00, m12 );
    Mout[1] = cuCadd( m01, m13 );
    Mout[2] = cuCadd( m20, m32 );
    Mout[3] = cuCadd( m21, m33 );
}

void conj2x2(cuDoubleComplex *M, cuDoubleComplex *Mout)
/* Calculate the conjugate of a matrix
 * It is safe for M and Mout to point to the same matrix
 */
{
    int i;
    for (i = 0; i < 4; i++)
        Mout[i] = cuConj(M[i]);
}


double norm2x2(cuDoubleComplex *M, cuDoubleComplex *Mout)
/* Normalise a 2x2 matrix via the Frobenius norm
 * It is safe for M and Mout to point to the same matrix.
 */
{
    // Calculate the normalising factor
    double Fnorm = 0.0;
    int i;
    for (i = 0; i < 4; i++)
        Fnorm += cuCreal( cuCmul( M[i], cuConj(M[i]) ) );

    Fnorm = sqrt(Fnorm);

    // Divide each element through by the normalising factor.
    // If norm is 0, then output zeros everywhere
    for (i = 0; i < 4; i++) {
        if (Fnorm == 0.0)
            Mout[i] = make_cuDoubleComplex( 0.0, 0.0 );
        else
            Mout[i] = make_cuDoubleComplex( cuCreal(M[i])/Fnorm, cuCimag(M[i])/Fnorm );
    }

    return Fnorm;
}


void dec2hms( char *out, double in, int sflag )
{
    int sign  = 1;
    char *ptr = out;
    int h, m;
    double s;

    if (in < 0.0)
    {
        sign = -1;
        in = fabs(in);
    }

    h = (int)in; in -= (double)h; in *= 60.0;
    m = (int)in; in -= (double)m; in *= 60.0;
    if (m >= 60)
    {
        // if minutes is 60 convert that to 1 hour
        h += 1;
        m -= 60;
    }
    s = in;
    if (s >= 59.995)
    {
        // if seconds is 60 convert that to 1 minute
        m += 1;
        s = 00.00;
    }
    if (sign==1 && sflag)
    {
        *ptr='+';
        ptr++;
    }
    else if (sign==-1)
    {
        *ptr='-';
        ptr++;
    }
    // Limiting the output's pointings' smallest significant figure to
    // 0.01 arc seconds
    sprintf( ptr, "%2.2d:%2.2d:%05.2f", h, m, s );
}



