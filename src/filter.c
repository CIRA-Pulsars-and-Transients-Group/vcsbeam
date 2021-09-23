/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cuComplex.h>

#include "vcsbeam.h"


cuDoubleComplex *roots_of_unity( int N )
/* Creates a complex-valued array containing the N roots of unity.
 * The caller must free this memory (by passing the returned pointer to
 * free()).
 */
{
    // Allocate memory
    cuDoubleComplex *roots = (cuDoubleComplex *)malloc( N*sizeof(cuDoubleComplex) );

    // Make sure memory was allocated correctly
    if (!roots)
    {
        fprintf( stderr, "error: roots_of_unity: could not allocate "
                         "memory\n" );
        exit(EXIT_FAILURE);
    }

    // Calculate the roots and store them in the roots array
    int n;
    for (n = 0; n < N; n++)
    {
        // e^{2πin/N}
        double th = 2*M_PI*(double)n/(double)N;
        roots[n] = make_cuDoubleComplex(cos(th), sin(th));
    }

    return roots;
}


pfb_filter *load_filter_coefficients( char *filtername, filter_type type, int nchans )
/* Load a set of filter coefficients
 * Inputs:
 *   FILTERNAME - string specifying a filter. There should be a corresponding
 *                file in the RUNTIME_DIR called FILTERNAME.dat.
 *   TYPE       - Whether it's an ANALYSIS_FILTER or a SYNTHESIS_FILTER
 *   NCHANS     - the number of channels that this filter will be applied to.
 *                For both ANALYSIS and SYNTHESIS filters, this should be
 *                the number of ANALYSIS channels.
 * Outputs:
 *   [return value] - pointer to the newly allocated struct containing the
 *                    coefficients. This should be freed using
 *                    free_pfb_filter()
 */
{
    // Make sure function arguments are not NULL
    if (filtername == NULL)
    {
        fprintf( stderr, "error: load_filter_coefficients: "
                "arguments must not be NULL\n" );
        exit(EXIT_FAILURE);
    }

    // Generate the file path from the given filtername
    char path[1024];
    sprintf( path, "%s/%s.dat", RUNTIME_DIR, filtername );

    // Open the file for reading
    FILE *f = fopen( path, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error: load_filter_coefficients: "
                "could not find filter '%s'\n", filtername );
        exit(EXIT_FAILURE);
    }

    // Go through file and count the number of coefficients
    // We assume that the file contains ONLY ASCII floating point
    // numbers separated by whitespace.
    double dummy;
    int num_read;

    pfb_filter *filter = (pfb_filter *)malloc( sizeof(pfb_filter) );
    filter->nchans     = nchans;
    filter->ncoeffs    = 0;
    filter->type       = type;

    while ((num_read = fscanf( f, "%lf", &dummy )) != EOF)
        filter->ncoeffs++;

    // Check to see that there is at least 1 coefficient!
    if (filter->ncoeffs == 0)
    {
        fprintf( stderr, "error: load_filter_coefficients: "
                "No coefficients found in '%s'\n", filtername );
        exit(EXIT_FAILURE);
    }
    filter->coeffs = (double *)malloc( filter->ncoeffs * sizeof(double) );

    // Rewind back to the beginning of the file and read them into
    // a new buffer.
    rewind( f );
    for (num_read = 0; num_read < filter->ncoeffs; num_read++)
        fscanf( f, "%lf", &(filter->coeffs[num_read]) );

    // Work out the number of taps, and issue a warning if the number
    // of channels does not divide evenly into the number of coefficients
    if (filter->ncoeffs % nchans != 0)
    {
        fprintf( stderr, "warning: load_filter_coefficients: "
                "number of channels (%d) does not divide evenly into the "
                "number of coefficients (%d)\n", nchans, filter->ncoeffs );
    }
    filter->ntaps = filter->ncoeffs / nchans;

    // Pre-calculate the twiddle factors
    filter->twiddles = roots_of_unity( nchans );

    // Close the file and return the pointer
    fclose( f );

    return filter;
}


void free_pfb_filter( pfb_filter *filter )
/* Free the memory allocated in load_filter_coefficients()
 */
{
    free( filter->coeffs );
    free( filter );
}