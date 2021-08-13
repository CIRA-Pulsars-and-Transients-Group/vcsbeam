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

void get_metafits_info( char *metafits, struct metafits_info *mi, unsigned int chan_width ) {
/* Read in the relevant information from the metafits file.
 * This function allocates dynamic memory. Destroy it with
 *   destroy_metafits_info(...)
 */

    fitsfile *fptr = NULL;
    int status     = 0;
    int anynull    = 0;

    // Open the metafits file
    fits_open_file(&fptr, metafits, READONLY, &status);
    if (fptr == NULL) {
        fprintf( stderr, "Failed to open metafits file \"%s\"\n", metafits );
        exit(EXIT_FAILURE);
    }

    // Read in the tile pointing information (and channel width)
    mi->date_obs = (char *)malloc( 32*sizeof(char) );
    fits_read_key(fptr, TDOUBLE, "RA",       &(mi->tile_pointing_ra),  NULL, &status);
    fits_read_key(fptr, TDOUBLE, "DEC",      &(mi->tile_pointing_dec), NULL, &status);
    fits_read_key(fptr, TDOUBLE, "AZIMUTH",  &(mi->tile_pointing_az),  NULL, &status);
    fits_read_key(fptr, TDOUBLE, "ALTITUDE", &(mi->tile_pointing_el),  NULL, &status);
    fits_read_key(fptr, TINT,    "EXPOSURE", &(mi->exposure),          NULL, &status);
    fits_read_key(fptr, TSTRING, "DATE-OBS", mi->date_obs,             NULL, &status);
    mi->chan_width = chan_width;

    if (status != 0) {
        fprintf(stderr, "Fits status set: failed to read az/alt, ");
        fprintf(stderr, "ra/dec fits keys from the header\n");
        exit(EXIT_FAILURE);
    }

    // Move to the binary table
    fits_movnam_hdu(fptr, BINARY_TBL, "TILEDATA", 0, &status);
    if (status != 0) {
        fprintf(stderr, "Error: Failed to move to TILEDATA HDU\n");
        exit(EXIT_FAILURE);
    }

    // Read in the number of inputs (= nstation * npol)
    fits_read_key(fptr, TINT, "NAXIS2", &(mi->ninput), NULL, &status);
    if (status != 0) {
        fprintf(stderr, "Error: Failed to read size of binary table in TILEDATA\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory
    mi->N_array         =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->E_array         =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->H_array         =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->cable_array     =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->flag_array      =       (int *)malloc( mi->ninput*sizeof(int)      );
    mi->weights_array   =    (double *)malloc( mi->ninput*sizeof(double)   );
    mi->antenna_num     = (short int *)malloc( mi->ninput*sizeof(short int));
    mi->tilenames       =     (char **)malloc( mi->ninput*sizeof(char *)   );
    mi->delays          =      (int **)malloc( mi->ninput*sizeof(int*)     );
    mi->amps            =   (double **)malloc( mi->ninput*sizeof(double*)  );
    int i;
    for (i = 0; i < (int)(mi->ninput); i++) {
        mi->tilenames[i] =   (char *)malloc( 32*sizeof(char)        );
        mi->delays[i]    =    (int *)malloc( NDELAYS*sizeof(int)    );
        mi->amps[i]      = (double *)malloc( NDELAYS*sizeof(double) );
    }
    char *testval = (char *) malloc(1024);

    /* Read the columns */
    int colnum;

    // Cable lengths
    for (i=0; i < (int)(mi->ninput); i++) {

        fits_get_colnum(fptr, 1, "Length", &colnum, &status);
        if(fits_read_col_str(fptr, colnum, i+1, 1, 1, "0.0", &testval, &anynull, &status)) {
            fprintf(stderr, "Error: Failed to cable column  in metafile\n");
            exit(EXIT_FAILURE);
        }

        sscanf(testval, "EL_%f", &(mi->cable_array[i]));
    }

    // Tile names
    fits_get_colnum(fptr, 1, "TileName", &colnum, &status);
    if (status != 0) {
        status = 0;
        fits_get_colnum(fptr, 1, "Tile", &colnum, &status);
    }
    if (status != 0) {
        fprintf(stderr, "Could not find either column \"TileName\" or \"Tile\" in metafits file\n");
        exit(EXIT_FAILURE);
    }

    fits_read_col(fptr, TSTRING, colnum, 1, 1, mi->ninput, NULL, mi->tilenames, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read Tile(Name) in metafile\n");
        exit(EXIT_FAILURE);
    }

    // North coordinate
    fits_get_colnum(fptr, 1, "North", &colnum, &status);
    fits_read_col_flt(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->N_array, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read  N coord in metafile\n");
        exit(EXIT_FAILURE);
    }

    // East coordinate
    fits_get_colnum(fptr, 1, "East", &colnum, &status);
    fits_read_col_flt(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->E_array, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read E coord in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Height coordinate
    fits_get_colnum(fptr, 1, "Height", &colnum, &status);
    fits_read_col_flt(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->H_array, &anynull, &status);

    if (status != 0){
        fprintf(stderr, "Error: Failed to read H coord in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Antenna number
    fits_get_colnum(fptr, 1, "Antenna", &colnum, &status);
    fits_read_col_sht(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->antenna_num, &anynull, &status);

    if (status != 0){
        fprintf(stderr, "Error: Failed to read field \"Antenna\" in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Flags & weights
    float wgt_sum = 0.0;
    fits_get_colnum(fptr, 1, "Flag", &colnum, &status);
    fits_read_col_int(fptr, colnum, 1, 1, mi->ninput, 0, mi->flag_array, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read flags column in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Beamforming delays (this only reads a single set of delays from the first tile,
    // and assumes that all tiles are the same)
    fits_get_colnum(fptr, 1, "Delays", &colnum, &status);
    for (i=0; i<mi->ninput; i++){
        fits_read_col_int(fptr, colnum, i+1, 1, NDELAYS, 0, mi->delays[i], &anynull, &status);
        if (status != 0){
            fprintf(stderr, "Error: Failed to read delays column in metafile\n");
            exit(EXIT_FAILURE);
        }

        // The amps should all be '1', except when the corresponding delay = '32'
        for (int j = 0; j < NDELAYS; j++){
            mi->amps[i][j] = (mi->delays[i][j] == 32 ? 0.0 : 1.0);
        }
    }
    // Invert value (flag off = full weight; flag on = zero weight)
    for (i = 0; i < mi->ninput; i++) {
        mi->weights_array[i] = 1.0 - (double)mi->flag_array[i];
        wgt_sum += mi->weights_array[i];
        // This differs from Ord's orig code, which sums squares. However,
        // all values should be = 1, so end result should be the same
    }

    // Exit with error if there are no weights
    if (wgt_sum == 0.0) {
        fprintf(stderr, "Zero weight sum on read\n");
        exit(EXIT_FAILURE);
    }

    // Clean up
    free( testval );
    fits_close_file(fptr, &status);
}


void destroy_metafits_info( struct metafits_info *mi ) {
/* Frees the memory allocated in the metafits_info struct
 */
    free( mi->N_array       );
    free( mi->E_array       );
    free( mi->H_array       );
    free( mi->cable_array   );
    free( mi->flag_array    );
    free( mi->weights_array );
    free( mi->antenna_num   );
    int i;
    for (i = 0; i < mi->ninput; i++){
        free( mi->tilenames[i] );
        free( mi->delays[i]    );
        free( mi->amps[i]      );
    }
    free( mi->tilenames     );
    free( mi->delays        );
    free( mi->amps          );
    free( mi->date_obs      );
}


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

void read_data( char *filename, uint8_t *data, int nbytes ) {

    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL) {
        fprintf( stderr, "Error opening file '%s'\n", filename );
        exit(EXIT_FAILURE);
    }

    // Read in nbytes bytes
    int nread = fread( (void *)data, sizeof(uint8_t), nbytes, f );

    // Check that I got all nbytes
    // Any other number is disallowed
    if (nread != nbytes) {
        fprintf( stderr, "Error: file %s does not contain %d bytes\n",
                filename, nbytes );
        exit(EXIT_FAILURE);
    }

    fclose(f);

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



