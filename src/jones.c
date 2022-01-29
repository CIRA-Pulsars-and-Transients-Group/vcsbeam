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

#include "vcsbeam.h"
#include "vcsbeam_private.h"

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



/*****************************
 * Generic matrix operations *
 *****************************/

/**
 * Copies a \f$2\times2\f$ complex-valued matrix.
 *
 * @param Min The source matrix, \f${\bf M}_\text{in}\f$.
 * @param[out] Mout The destination matrix, \f${\bf M}_\text{out} = {\bf M}_\text{in}\f$.
 */
void cp2x2(cuDoubleComplex *Min, cuDoubleComplex *Mout)
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
cuDoubleComplex reciprocal_complex( cuDoubleComplex z )
{
    double scale = 1.0/(z.x*z.x + z.y*z.y);
    return make_cuDoubleComplex( scale*z.x, -scale*z.y );
}

/**
 * Computes the negative of a complex number.
 *
 * @param z A complex number
 * @return \f$-z\f$
 */
cuDoubleComplex negate_complex( cuDoubleComplex z )
{
    return make_cuDoubleComplex( -z.x, -z.y );
}

/**
 * Computes the inverse of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param Min The input matrix, \f${\bf M}_\text{in}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}_\text{in}^{-1}\f$
 */
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


/**
 * Performs a matrix multiplication of two \f$2\times2\f$ complex-valued matrices.
 *
 * @param M1 The first input matrix, \f${\bf M}_1\f$
 * @param M2 The second input matrix, \f${\bf M}_2\f$
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}_1 \times {\bf M}_2\f$.
 */
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

/**
 * Performs a matrix multiplication of a \f$2\times2\f$ real-valued matrix
 * (on the left) with a \f$2\times2\f$ complex-valued matrix (on the right).
 *
 * @param M1 The first input matrix (real-valued), \f${\bf M}_1\f$
 * @param M2 The second input matrix (complex-valued), \f${\bf M}_2\f$
 * @param[out] Mout The output matrix (complex-valued), \f${\bf M}_\text{out} = {\bf M}_1 \times {\bf M}_2\f$.
 */
void mult2x2d_RxC(double *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout)
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

/**
 * Performs a matrix multiplication of a \f$2\times2\f$ complex-valued matrix
 * (on the left) with a \f$2\times2\f$ real-valued matrix (on the right).
 *
 * @param M1 The first input matrix (complex-valued), \f${\bf M}_1\f$
 * @param M2 The second input matrix (real-valued), \f${\bf M}_2\f$
 * @param[out] Mout The output matrix (complex-valued), \f${\bf M}_\text{out} = {\bf M}_1 \times {\bf M}_2\f$.
 */
void mult2x2d_CxR( cuDoubleComplex *M1, double *M2, cuDoubleComplex *Mout )
{
    cuDoubleComplex m00 = make_cuDoubleComplex( cuCreal(M1[0])*M2[0], cuCimag(M1[0])*M2[0] );
    cuDoubleComplex m12 = make_cuDoubleComplex( cuCreal(M1[1])*M2[2], cuCimag(M1[1])*M2[2] );
    cuDoubleComplex m01 = make_cuDoubleComplex( cuCreal(M1[0])*M2[1], cuCimag(M1[0])*M2[1] );
    cuDoubleComplex m13 = make_cuDoubleComplex( cuCreal(M1[1])*M2[3], cuCimag(M1[1])*M2[3] );
    cuDoubleComplex m20 = make_cuDoubleComplex( cuCreal(M1[2])*M2[0], cuCimag(M1[2])*M2[0] );
    cuDoubleComplex m32 = make_cuDoubleComplex( cuCreal(M1[3])*M2[2], cuCimag(M1[3])*M2[2] );
    cuDoubleComplex m21 = make_cuDoubleComplex( cuCreal(M1[2])*M2[1], cuCimag(M1[2])*M2[1] );
    cuDoubleComplex m33 = make_cuDoubleComplex( cuCreal(M1[3])*M2[3], cuCimag(M1[3])*M2[3] );
    Mout[0] = cuCadd( m00, m12 );
    Mout[1] = cuCadd( m01, m13 );
    Mout[2] = cuCadd( m20, m32 );
    Mout[3] = cuCadd( m21, m33 );
}

/**
 * Computes the conjugate of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = {\bf M}^\ast\f$.
 *
 * It is safe for `M` and `Mout` to point to the same matrix.
 */
void conj2x2(cuDoubleComplex *M, cuDoubleComplex *Mout)
{
    int i;
    for (i = 0; i < 4; i++)
        Mout[i] = cuConj(M[i]);
}


/**
 * Normalises a \f$2\times2\f$ matrix via the Frobenius norm.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] Mout The output matrix, \f${\bf M}_\text{out} = |{\bf M}_\text{in}|\f$.
 *
 * It is safe for `M` and `Mout` to point to the same matrix.
 */
double norm2x2(cuDoubleComplex *M, cuDoubleComplex *Mout)
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
void reverse2x2( cuDoubleComplex *M, cuDoubleComplex *Mout )
{
    cuDoubleComplex m0 = M[0];
    cuDoubleComplex m1 = M[1];
    cuDoubleComplex m2 = M[2];
    cuDoubleComplex m3 = M[3];

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
void swaprows2x2( cuDoubleComplex *M, cuDoubleComplex *Mout )
{
    cuDoubleComplex m0 = M[0];
    cuDoubleComplex m1 = M[1];
    cuDoubleComplex m2 = M[2];
    cuDoubleComplex m3 = M[3];

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
void swapcols2x2( cuDoubleComplex *M, cuDoubleComplex *Mout )
{
    cuDoubleComplex m0 = M[0];
    cuDoubleComplex m1 = M[1];
    cuDoubleComplex m2 = M[2];
    cuDoubleComplex m3 = M[3];

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
bool is2x2zero( cuDoubleComplex *M )
{
    return (cuCreal(M[0]) == 0.0 &&
            cuCimag(M[0]) == 0.0 &&
            cuCreal(M[1]) == 0.0 &&
            cuCimag(M[1]) == 0.0 &&
            cuCreal(M[2]) == 0.0 &&
            cuCimag(M[2]) == 0.0 &&
            cuCreal(M[3]) == 0.0 &&
            cuCimag(M[3]) == 0.0);
}

/**
 * Computes the Hermitian (i.e. conjugate transpose) of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] H The output matrix, \f${\bf H} = {\bf M}^\dagger\f$.
 *
 * It is safe for `M` and `H` to point to the same matrix.
 */
void calc_hermitian( cuDoubleComplex *M, cuDoubleComplex *H )
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

/**
 * Computes the coherency matrix of a \f$2\times2\f$ complex-valued matrix.
 *
 * @param M The input matrix, \f${\bf M}\f$.
 * @param[out] C The output matrix, \f${\bf C} = {\bf M}{\bf M}^\dagger\f$.
 *
 * It is safe for `M` and `C` to point to the same matrix.
 */
void calc_coherency_matrix( cuDoubleComplex *M, cuDoubleComplex *C )
/* Compute C = M x M^H, where "^H" is the hermitian operation
 */
{
    cuDoubleComplex MH[4];
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
void fprintf_complex_matrix( FILE *fout, cuDoubleComplex *M )
{
    fprintf( fout, "[ %lf%+lf*i, %lf%+lf*i; %lf%+lf*i, %lf%+lf*i ]\n",
                   cuCreal(M[0]), cuCimag(M[0]),
                   cuCreal(M[1]), cuCimag(M[1]),
                   cuCreal(M[2]), cuCimag(M[2]),
                   cuCreal(M[3]), cuCimag(M[3])
           );
}
