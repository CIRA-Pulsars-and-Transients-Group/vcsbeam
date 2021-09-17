/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <cuComplex.h>

#include <mwalib.h>
#include <star/pal.h>
#include <star/palmac.h>
#include <psrfits.h>
#include <fitsio.h>
#include <mwa_hyperbeam.h>

#include "vcsbeam.h"
#include "primary_beam.h"

void create_antenna_lists( MetafitsMetadata *obs_metadata, uint32_t *polX_idxs, uint32_t *polY_idxs )
/* Creates a list of indexes into the data for the X and Y polarisations,
 * ordered by antenna number. Assumption: polX_idxs and polY_idxs have
 * sufficient allocated memory.
 */
{
    // Go through the rf_inputs and construct the lookup table for the antennas
    unsigned int ninputs = obs_metadata->num_rf_inputs;
    unsigned int i, ant;
    for (i = 0; i < ninputs; i++)
    {
        ant = obs_metadata->rf_inputs[i].ant;
        if (*(obs_metadata->rf_inputs[i].pol) == 'X')
            polX_idxs[ant] = i;
        else // if (*(obs_metadata->rf_inputs.pol) == 'Y')
            polY_idxs[ant] = i;
    }
}

void get_jones(
        // an array of pointings [pointing][ra/dec][characters]
        int                    npointing, // number of pointings
        MetafitsMetadata      *obs_metadata,
        int                    coarse_chan_idx,
        struct                 calibration *cal,
        cuDoubleComplex       *D,
        cuDoubleComplex       *B,
        cuDoubleComplex       *invJi )
{

    // Give "shorthand" variables for often-used values in metafits
    long int frequency = obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_start_hz;
    int nant           = obs_metadata->num_ants;
    int nchan          = obs_metadata->num_volt_fine_chans_per_coarse;
    int npol           = obs_metadata->num_ant_pols;   // (X,Y)
    int chan_width     = obs_metadata->corr_fine_chan_width_hz;

    int p;       // Pointing number
    int ant;     // Antenna number
    int ch;      // Channel number
    int p1, p2;  // Counters for polarisation

    /* easy -- now the positions from the database */

    cuDoubleComplex Ji[npol*npol];              // Gain in Desired Direction

    long int freq_ch;

    double uv_angle;
    cuDoubleComplex uv_phase; // For the UV phase correction

    double Fnorm;

    int j_idx;

    for (p = 0; p < npointing; p++)
    {
        // Everything from this point on is frequency-dependent
        for (ch = 0; ch < nchan; ch++) {

            // Calculating direction-dependent matrices
            freq_ch = frequency + ch*chan_width;    // The frequency of this fine channel

            // Calculate the UV phase correction for this channel
            uv_angle = cal->phase_slope*freq_ch + cal->phase_offset;
            uv_phase = make_cuDoubleComplex( cos(uv_angle), sin(uv_angle) );

            for (ant = 0; ant < nant; ant++)
            {
                // The index to the first element in the Jones matrix for this
                // antenna and channel. Applies to both the D and J arrays.
                j_idx = J_IDX(ant,ch,0,0,nchan,npol);

                mult2x2d(&(D[j_idx]), &(B[PB_IDX(p, ant, 0, nant, npol*npol)]), Ji); // the gain in the desired look direction

                // Apply the UV phase correction (to the bottom row of the Jones matrix)
                Ji[2] = cuCmul( Ji[2], uv_phase );
                Ji[3] = cuCmul( Ji[3], uv_phase );

                // Now, calculate the inverse Jones matrix

                // Ord's original comment for the following line is:
                // "The RTS conjugates the sky so beware"
                conj2x2( Ji, Ji );

                Fnorm = norm2x2( Ji, Ji );

                if (Fnorm != 0.0)
                    inv2x2S( Ji, &(invJi[j_idx]) );
                else {
                    for (p1 = 0; p1 < npol;  p1++)
                    for (p2 = 0; p2 < npol;  p2++)
                        invJi[J_IDX(ant,ch,p1,p2,nchan,npol)] = make_cuDoubleComplex( 0.0, 0.0 );
                }

            } // end loop through antenna/pol (rf_input)
        } // end loop through fine channels (ch)
    } // end loop through pointings (p)
}


/*****************************
 * Generic matrix operations *
 *****************************/

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


void inv2x2S(cuDoubleComplex *Min, cuDoubleComplex *Mout)
// Same as inv2x2(), but the output is a 2x2 2D array, instead of a 4-element
// 1D array
{
    cuDoubleComplex m1 = cuCmul( Min[0], Min[3] );
    cuDoubleComplex m2 = cuCmul( Min[1], Min[2] );
    cuDoubleComplex det = cuCsub( m1, m2 );
    cuDoubleComplex inv_det = reciprocal_complex( det );
    Mout[0] = cuCmul(       inv_det,  Min[3] );
    Mout[1] = cuCmul( negate_complex(inv_det), Min[1] );
    Mout[2] = cuCmul( negate_complex(inv_det), Min[2] );
    Mout[3] = cuCmul(       inv_det,  Min[0] );
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


void calc_hermitian( cuDoubleComplex *M, cuDoubleComplex *H )
/* Calculate H = M^H, where "^H" is the hermitian operation
 */
{
    // Put in temporary matrix, so that H may point to the same memory
    // address as M
    cuDoubleComplex tmp[4];
    tmp[0] = cuConj( M[0] );
    tmp[1] = cuConj( M[2] );
    tmp[2] = cuConj( M[1] );
    tmp[3] = cuConj( M[3] );

    cp2x2( tmp, H );
}

void calc_coherency_matrix( cuDoubleComplex *M, cuDoubleComplex *C )
/* Compute C = M x M^H, where "^H" is the hermitian operation
 */
{
    cuDoubleComplex MH[4];
    calc_hermitian( M, MH );
    mult2x2d( M, MH, C );
}

