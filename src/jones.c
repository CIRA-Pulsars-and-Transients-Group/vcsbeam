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

#include <mwalib.h>

#include "vcsbeam.h"
#include "gpu_macros.h"

/**
 * Creates arrays of indexes for antennas and polarisation according to the
 * ordering used in legacy VCS observations.
 *
 * @param vm The VCSBeam context struct
 *
 * The index arrays are needed for pairing the antenna metadata with the
 * voltage samples, which historically were ordered in a different way.
 * Fortunately, the metafits file contains the VCS ordering in the "VCSOrder"
 * field. This function converts these indexes into a lookup array that can
 * be conveniently passed to GPU kernels at runtime. A separate array is
 * made for each polarisation, with the indexes into the arrays being the
 * metafits "Antenna" number, and the value stored in that array position
 * being the "VCSOrder".
 *
 * The metafits information is drawn from `vm&rarr;obs_metadata`, and the
 * indexes are saved in `vm&rarr;polP_idxs` and `vm&rarr;polQ_idxs`.
 *
 * No memory is allocated by this function. It is assumed that all arrays
 * have been previously allocated, and are sufficiently large.
 */
void vmSetPolIdxLists( vcsbeam_context *vm )
{
    // Go through the rf_inputs and construct the lookup table for the antennas
    unsigned int ninputs = vm->obs_metadata->num_rf_inputs;
    unsigned int ant;
    unsigned int i;
    char pol;
    for (i = 0; i < ninputs; i++)
    {
        ant = vm->obs_metadata->rf_inputs[i].ant;
        pol = *(vm->obs_metadata->rf_inputs[i].pol);

        // As written in the documentation, the polarisation identified by
        // 'X' in the metafits file is the N-S-aligned, or 'Q' dipole,
        // 'Y' in the metafits file is the E-W-aligned, or 'P' dipole,
        if (pol == 'X')
        {
            vm->polQ_idxs[ant] = vm->obs_metadata->rf_inputs[i].vcs_order;
        }
        else // if (pol == 'Y')
        {
            vm->polP_idxs[ant] = vm->obs_metadata->rf_inputs[i].vcs_order;
        }
    }
}

/**
 * Computes the (inverse) Jones matrix \f${\bf J}^{-1} = \left({\bf D}{\bf B}\right)^{-1}\f$.
 *
 * \f${\bf J}^{-1}\f$ is the Jones matrix that converts from instrumental
 * voltages to sky voltages, according to
 * \f[
 * \begin{aligned}
 *     {\bf e} &= {\bf J}^{-1} {\bf v}, \\
 *     \begin{bmatrix} e_x \\ e_y \end{bmatrix} &=
 *     \begin{bmatrix} J_{xq} & J_{xp} \\ J_{yq} & J_{yp} \end{bmatrix}
 *     \begin{bmatrix} v_q \\ v_p \end{bmatrix}.
 * \end{aligned}
 * \f]
 * It is computed by multiplying the instrumental calibration Jones matrix
 * with the beam model matrix, and taking the inverse
 * \f[{\bf J}^{-1} = \left({\bf D}{\bf B}\right)^{-1},\f]
 * where \f${\bf D}\f$ is taken from `vm&rarr;D`, and \f${\bf B}\f$ from
 * `vm&rarr;pb.B`.
 * The result is stored in `vm&rarr;J`.
 *
 * \f${\bf J}^{-1}\f$ has dimensions \f$(N_a \times N_f \times N_p \times N_p)\f$,
 * where
 *  - \f$N_a\f$ is the number of antennas (`vm&rarr;obs_metadata&rarr;num_ants`)
 *  - \f$N_f\f$ is the number of frequencies (`vm&rarr;vm&rarr;nfine_chan`)
 *  - \f$N_p\f$ is the number of polarisation (`vm&rarr;obs_metadata&rarr;num_ant_pols`)
 *
 * The antenna ordering is given by the "Antenna" keyword in the observation's
 * metafits file. The frequencies are ordered from lowest to highest.
 *
 * @todo <a href="https://github.com/CIRA-Pulsars-and-Transients-Group/vcsbeam/issues/9">Issue #9</a>
 */
void vmCalcJ( vcsbeam_context *vm )
{

    // Give "shorthand" variables for often-used values in metafits
    int nant           = vm->obs_metadata->num_ants;
    int nchan          = vm->nfine_chan;
    int npol           = vm->obs_metadata->num_ant_pols;   // (X,Y)

    unsigned int p;  // Pointing number
    int ant;         // Antenna number
    int ch;          // Channel number
    int p1, p2;      // Counters for polarisation

    /* easy -- now the positions from the database */

    gpuDoubleComplex Ji[npol*npol];              // Gain in Desired Direction

    double Fnorm;

    int d_idx, j_idx, pb_idx;

    for (p = 0; p < vm->npointing; p++)
    {
        // Everything from this point on is frequency-dependent
        for (ch = 0; ch < nchan; ch++) {

            for (ant = 0; ant < nant; ant++)
            {
                // The index to the first element in the Jones matrix for this
                // antenna and channel. Applies to both the D and J arrays.
                d_idx  = D_IDX(ant,ch,0,0,nchan,npol);
                j_idx  = J_IDX(p,ant,ch,0,0,nant,nchan,npol);
                pb_idx = PB_IDX(p, ant, 0, nant, npol*npol);

                mult2x2d(&(vm->D[d_idx]), &(vm->pb.B[pb_idx]), Ji); // the gain in the desired look direction

#ifdef DEBUG
if (ch == 50 && ant == 0)
{
    fprintf( stderr, "pointing = %u:\n", p );
    fprintf( stderr, "\tD       = "); fprintf_complex_matrix( stderr, &(vm->D[d_idx]) );
    fprintf( stderr, "\tBP      = "); fprintf_complex_matrix( stderr, &(vm->pb.B[pb_idx]) );
    fprintf( stderr, "\tJ = DBP = "); fprintf_complex_matrix( stderr, Ji );
}
#endif

                // Now, calculate the inverse Jones matrix
                Fnorm = norm2x2( Ji, Ji );

                if (Fnorm != 0.0)
                    inv2x2S( Ji, &(vm->J[j_idx]) );
                else {
                    for (p1 = 0; p1 < npol;  p1++)
                    for (p2 = 0; p2 < npol;  p2++)
                        vm->J[J_IDX(p,ant,ch,p1,p2,nant,nchan,npol)] = make_gpuDoubleComplex( 0.0, 0.0 );
                }

            } // end loop through antenna/pol (rf_input)
        } // end loop through fine channels (ch)
    } // end loop through pointings (p)
}

/**
 * Wrapper function for vmCalcPhi(), vmCalcB(), and vmCalcJ().
 *
 * @param vm The VCSBeam context struct
 * @param ras_hours An array of right ascensions, in decimal hours, one for each pointing.
 * @param decs_degs An array of declinations, in decimal degrees, one for each pointing.
 * @param beam_geom_vals An array of beam geometries, one for each pointing.
 *
 * For each pointing, the quantities \f$e^{i\varphi}\f$, \f${\bf B}\f$, and \f${\bf J}^{-1}\f$ are calculated.
 * (\f${\bf D}\f$ is not recalculated, but is used in the calculation of \f${\bf J}^{-1}\f$).
 */
void vmCalcJonesAndDelays( vcsbeam_context *vm, double *ras_hours, double *decs_degs, beam_geom *beam_geom_vals )
{
    logger_start_stopwatch( vm->log, "delay", true );

    uintptr_t timestep_idx = vm->chunk_to_load / vm->chunks_per_second;
    double sec_offset = (double)(timestep_idx + vm->gps_seconds_to_process[0] - vm->obs_metadata->obs_id);
    double mjd = vm->obs_metadata->sched_start_mjd + (sec_offset + 0.5)/86400.0;

    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
        calc_beam_geom( ras_hours[p], decs_degs[p], mjd, &beam_geom_vals[p] );

    // Calculate the geometric delays, primary beam and Jones matrices
    vmCalcPhi( vm, beam_geom_vals );
    vmCalcB( vm, beam_geom_vals );
    vmCalcJ( vm );

    logger_stop_stopwatch( vm->log, "delay" );
}


/*****************************
 * Generic matrix operations *
 *****************************/

/**
 * Copies a \f$2\times2\f$ complex-valued matrix.
 *
 * @param Min The source matrix, \f${\bf M}_\text{in}\f$.
 * @param[out] Mout The destination matrix, \f${\bf M}_\text{out} = {\bf M}_\text{in}\f$.
 */
void cp2x2(gpuDoubleComplex *Min, gpuDoubleComplex *Mout)
{
    Mout[0] = Min[0];
    Mout[1] = Min[1];
    Mout[2] = Min[2];
    Mout[3] = Min[3];
}

/**
 * Computes the reciprocal of a complex number.
 *
 * @param z A complex number
 * @return \f$\dfrac{1}{z}\f$
 */
gpuDoubleComplex reciprocal_complex( gpuDoubleComplex z )
{
    double scale = 1.0/(z.x*z.x + z.y*z.y);
    return make_gpuDoubleComplex( scale*z.x, -scale*z.y );
}

/**
 * Computes the negative of a complex number.
 *
 * @param z A complex number
 * @return \f$-z\f$
 */
gpuDoubleComplex negate_complex( gpuDoubleComplex z )
{
    return make_gpuDoubleComplex( -z.x, -z.y );
}

/**
 * Computes the inverse of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param Min The input matrix, \f${\bf M}_\text{in}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}_\text{in}^{-1}\f$
 */
void inv2x2(gpuDoubleComplex *Min, gpuDoubleComplex *Mout)
{
    gpuDoubleComplex m00 = Min[0];
    gpuDoubleComplex m01 = Min[1];
    gpuDoubleComplex m10 = Min[2];
    gpuDoubleComplex m11 = Min[3];

    gpuDoubleComplex m1 = gpuCmul( m00, m11 );
    gpuDoubleComplex m2 = gpuCmul( m01, m10 );

    gpuDoubleComplex det = gpuCsub( m1, m2 );
    gpuDoubleComplex inv_det = reciprocal_complex( det );

    Mout[0] = gpuCmul(       inv_det,  m11 );
    Mout[1] = gpuCmul( negate_complex(inv_det), m01 );
    Mout[2] = gpuCmul( negate_complex(inv_det), m10 );
    Mout[3] = gpuCmul(       inv_det,  m00 );
}

/**
 * Computes the inverse of a \f$2\times2\f$ real-valued matrix.
 *
 * @param Min The input matrix, \f${\bf M}_\text{in}\f$
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}_\text{in}^{-1}\f$.
 */
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


/**
 * Computes the inverse of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param Min The input matrix, \f${\bf M}_\text{in}\f$
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}_\text{in}^{-1}\f$
 *
 * \see inv2x2()
 *
 * @todo Does this do the same thing as inv2x2()? Check the code, and delete
 *       if not needed.
 */
void inv2x2S(gpuDoubleComplex *Min, gpuDoubleComplex *Mout)
// Same as inv2x2(), but the output is a 2x2 2D array, instead of a 4-element
// 1D array
{
    gpuDoubleComplex m1 = gpuCmul( Min[0], Min[3] );
    gpuDoubleComplex m2 = gpuCmul( Min[1], Min[2] );
    gpuDoubleComplex det = gpuCsub( m1, m2 );
    gpuDoubleComplex inv_det = reciprocal_complex( det );
    Mout[0] = gpuCmul(       inv_det,  Min[3] );
    Mout[1] = gpuCmul( negate_complex(inv_det), Min[1] );
    Mout[2] = gpuCmul( negate_complex(inv_det), Min[2] );
    Mout[3] = gpuCmul(       inv_det,  Min[0] );
}


/**
 * Performs a matrix multiplication of two \f$2\times2\f$ complex-valued matrices.
 *
 * @param M1 The first input matrix, \f${\bf M}_1\f$
 * @param M2 The second input matrix, \f${\bf M}_2\f$
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}_1 \times {\bf M}_2\f$.
 */
void mult2x2d(gpuDoubleComplex *M1, gpuDoubleComplex *M2, gpuDoubleComplex *Mout)
{
    gpuDoubleComplex m00 = gpuCmul( M1[0], M2[0] );
    gpuDoubleComplex m12 = gpuCmul( M1[1], M2[2] );
    gpuDoubleComplex m01 = gpuCmul( M1[0], M2[1] );
    gpuDoubleComplex m13 = gpuCmul( M1[1], M2[3] );
    gpuDoubleComplex m20 = gpuCmul( M1[2], M2[0] );
    gpuDoubleComplex m32 = gpuCmul( M1[3], M2[2] );
    gpuDoubleComplex m21 = gpuCmul( M1[2], M2[1] );
    gpuDoubleComplex m33 = gpuCmul( M1[3], M2[3] );
    Mout[0] = gpuCadd( m00, m12 );
    Mout[1] = gpuCadd( m01, m13 );
    Mout[2] = gpuCadd( m20, m32 );
    Mout[3] = gpuCadd( m21, m33 );
}

/**
 * Performs a matrix multiplication of a \f$2\times2\f$ real-valued matrix
 * (on the left) with a \f$2\times2\f$ complex-valued matrix (on the right).
 *
 * @param M1 The first input matrix (real-valued), \f${\bf M}_1\f$
 * @param M2 The second input matrix (complex-valued), \f${\bf M}_2\f$
 * @param[out] Mout The output matrix (complex-valued), \f${\bf M}_\text{out} = {\bf M}_1 \times {\bf M}_2\f$.
 */
void mult2x2d_RxC(double *M1, gpuDoubleComplex *M2, gpuDoubleComplex *Mout)
{
    gpuDoubleComplex m00 = make_gpuDoubleComplex( M1[0]*gpuCreal(M2[0]), M1[0]*gpuCimag(M2[0]) );
    gpuDoubleComplex m12 = make_gpuDoubleComplex( M1[1]*gpuCreal(M2[2]), M1[1]*gpuCimag(M2[2]) );
    gpuDoubleComplex m01 = make_gpuDoubleComplex( M1[0]*gpuCreal(M2[1]), M1[0]*gpuCimag(M2[1]) );
    gpuDoubleComplex m13 = make_gpuDoubleComplex( M1[1]*gpuCreal(M2[3]), M1[1]*gpuCimag(M2[3]) );
    gpuDoubleComplex m20 = make_gpuDoubleComplex( M1[2]*gpuCreal(M2[0]), M1[2]*gpuCimag(M2[0]) );
    gpuDoubleComplex m32 = make_gpuDoubleComplex( M1[3]*gpuCreal(M2[2]), M1[3]*gpuCimag(M2[2]) );
    gpuDoubleComplex m21 = make_gpuDoubleComplex( M1[2]*gpuCreal(M2[1]), M1[2]*gpuCimag(M2[1]) );
    gpuDoubleComplex m33 = make_gpuDoubleComplex( M1[3]*gpuCreal(M2[3]), M1[3]*gpuCimag(M2[3]) );
    Mout[0] = gpuCadd( m00, m12 );
    Mout[1] = gpuCadd( m01, m13 );
    Mout[2] = gpuCadd( m20, m32 );
    Mout[3] = gpuCadd( m21, m33 );
}

/**
 * Performs a matrix multiplication of a \f$2\times2\f$ complex-valued matrix
 * (on the left) with a \f$2\times2\f$ real-valued matrix (on the right).
 *
 * @param M1 The first input matrix (complex-valued), \f${\bf M}_1\f$
 * @param M2 The second input matrix (real-valued), \f${\bf M}_2\f$
 * @param[out] Mout The output matrix (complex-valued), \f${\bf M}_\text{out} = {\bf M}_1 \times {\bf M}_2\f$.
 */
void mult2x2d_CxR( gpuDoubleComplex *M1, double *M2, gpuDoubleComplex *Mout )
{
    gpuDoubleComplex m00 = make_gpuDoubleComplex( gpuCreal(M1[0])*M2[0], gpuCimag(M1[0])*M2[0] );
    gpuDoubleComplex m12 = make_gpuDoubleComplex( gpuCreal(M1[1])*M2[2], gpuCimag(M1[1])*M2[2] );
    gpuDoubleComplex m01 = make_gpuDoubleComplex( gpuCreal(M1[0])*M2[1], gpuCimag(M1[0])*M2[1] );
    gpuDoubleComplex m13 = make_gpuDoubleComplex( gpuCreal(M1[1])*M2[3], gpuCimag(M1[1])*M2[3] );
    gpuDoubleComplex m20 = make_gpuDoubleComplex( gpuCreal(M1[2])*M2[0], gpuCimag(M1[2])*M2[0] );
    gpuDoubleComplex m32 = make_gpuDoubleComplex( gpuCreal(M1[3])*M2[2], gpuCimag(M1[3])*M2[2] );
    gpuDoubleComplex m21 = make_gpuDoubleComplex( gpuCreal(M1[2])*M2[1], gpuCimag(M1[2])*M2[1] );
    gpuDoubleComplex m33 = make_gpuDoubleComplex( gpuCreal(M1[3])*M2[3], gpuCimag(M1[3])*M2[3] );
    Mout[0] = gpuCadd( m00, m12 );
    Mout[1] = gpuCadd( m01, m13 );
    Mout[2] = gpuCadd( m20, m32 );
    Mout[3] = gpuCadd( m21, m33 );
}

/**
 * Computes the conjugate of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}^\ast\f$.
 *
 * It is safe for `M` and `Mout` to point to the same matrix.
 */
void conj2x2(gpuDoubleComplex *M, gpuDoubleComplex *Mout)
{
    int i;
    for (i = 0; i < 4; i++)
        Mout[i] = gpuConj(M[i]);
}


/**
 * Normalises a \f$2\times2\f$ matrix via the Frobenius norm.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = |{\bf M}_\text{in}|\f$.
 *
 * It is safe for `M` and `Mout` to point to the same matrix.
 */
double norm2x2(gpuDoubleComplex *M, gpuDoubleComplex *Mout)
{
    // Calculate the normalising factor
    double Fnorm = 0.0;
    int i;
    for (i = 0; i < 4; i++)
        Fnorm += gpuCreal( gpuCmul( M[i], gpuConj(M[i]) ) );

    Fnorm = sqrt(Fnorm);

    // Divide each element through by the normalising factor.
    // If norm is 0, then output zeros everywhere
    for (i = 0; i < 4; i++) {
        if (Fnorm == 0.0)
            Mout[i] = make_gpuDoubleComplex( 0.0, 0.0 );
        else
            Mout[i] = make_gpuDoubleComplex( gpuCreal(M[i])/Fnorm, gpuCimag(M[i])/Fnorm );
    }

    return Fnorm;
}


/**
 * Reverses the elements of a \f$2\times2\f$ matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out}\f$.
 *
 * \f[\begin{bmatrix} 0 & 1 \\ 2 & 3\end{bmatrix}
 * \rightarrow
 * \begin{bmatrix} 3 & 2 \\ 1 & 0 \end{bmatrix}\f]
 *
 * It is safe for `M` and `Mout` to point to the same matrix.
 */
void reverse2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout )
{
    gpuDoubleComplex m0 = M[0];
    gpuDoubleComplex m1 = M[1];
    gpuDoubleComplex m2 = M[2];
    gpuDoubleComplex m3 = M[3];

    Mout[0] = m3;
    Mout[1] = m2;
    Mout[2] = m1;
    Mout[3] = m0;
}

/**
 * Swaps the rows of a \f$2\times2\f$ matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out}\f$.
 *
 * \f[\begin{bmatrix} 0 & 1 \\ 2 & 3\end{bmatrix}
 * \rightarrow
 * \begin{bmatrix} 2 & 3 \\ 0 & 1 \end{bmatrix}\f]
 *
 * It is safe for `M` and `Mout` to point to the same matrix.
 */
void swaprows2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout )
{
    gpuDoubleComplex m0 = M[0];
    gpuDoubleComplex m1 = M[1];
    gpuDoubleComplex m2 = M[2];
    gpuDoubleComplex m3 = M[3];

    Mout[0] = m2;
    Mout[1] = m3;
    Mout[2] = m0;
    Mout[3] = m1;
}

/**
 * Swaps the columns of a \f$2\times2\f$ matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out}\f$.
 *
 * \f[\begin{bmatrix} 0 & 1 \\ 2 & 3\end{bmatrix}
 * \rightarrow
 * \begin{bmatrix} 1 & 0 \\ 3 & 2 \end{bmatrix}\f]
 *
 * It is safe for `M` and `Mout` to point to the same matrix.
 */
void swapcols2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout )
{
    gpuDoubleComplex m0 = M[0];
    gpuDoubleComplex m1 = M[1];
    gpuDoubleComplex m2 = M[2];
    gpuDoubleComplex m3 = M[3];

    Mout[0] = m1;
    Mout[1] = m0;
    Mout[2] = m3;
    Mout[3] = m2;
}

/**
 * Returns true if and only if all elements of a \f$2\times2\f$ complex-valued matrix are identically zero.
 *
 * @param M The matrix to be tested, \f${\bf M}\f$.
 */
bool is2x2zero( gpuDoubleComplex *M )
{
    return (gpuCreal(M[0]) == 0.0 &&
            gpuCimag(M[0]) == 0.0 &&
            gpuCreal(M[1]) == 0.0 &&
            gpuCimag(M[1]) == 0.0 &&
            gpuCreal(M[2]) == 0.0 &&
            gpuCimag(M[2]) == 0.0 &&
            gpuCreal(M[3]) == 0.0 &&
            gpuCimag(M[3]) == 0.0);
}

/**
 * Computes the Hermitian (i.e. conjugate transpose) of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] H The output matrix, \f${\bf H} = {\bf M}^\dagger\f$.
 *
 * It is safe for `M` and `H` to point to the same matrix.
 */
void calc_hermitian( gpuDoubleComplex *M, gpuDoubleComplex *H )
{
    // Put in temporary matrix, so that H may point to the same memory
    // address as M
    gpuDoubleComplex tmp[4];
    tmp[0] = gpuConj( M[0] );
    tmp[1] = gpuConj( M[2] );
    tmp[2] = gpuConj( M[1] );
    tmp[3] = gpuConj( M[3] );

    cp2x2( tmp, H );
}

/**
 * Computes the coherency matrix of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] C The output matrix, \f${\bf C} = {\bf M}{\bf M}^\dagger\f$.
 *
 * It is safe for `M` and `C` to point to the same matrix.
 */
void calc_coherency_matrix( gpuDoubleComplex *M, gpuDoubleComplex *C )
/* Compute C = M x M^H, where "^H" is the hermitian operation
 */
{
    gpuDoubleComplex MH[4];
    calc_hermitian( M, MH );
    mult2x2d( M, MH, C );
}


/**
 * Prints a \f$2\times2\f$ complex-valued matrix.
 *
 * @param fout The file stream to output to
 * @param M The matrix to print
 *
 * The output format is
 *
 * `[ a+b*i, c+d*i; e+f*i, g+h*i ]`
 *
 * and a newline is printed at the end.
 */
void fprintf_complex_matrix( FILE *fout, gpuDoubleComplex *M )
{
    fprintf( fout, "[ %lf%+lf*i, %lf%+lf*i; %lf%+lf*i, %lf%+lf*i ]\n",
                   gpuCreal(M[0]), gpuCimag(M[0]),
                   gpuCreal(M[1]), gpuCimag(M[1]),
                   gpuCreal(M[2]), gpuCimag(M[2]),
                   gpuCreal(M[3]), gpuCimag(M[3])
           );
}
