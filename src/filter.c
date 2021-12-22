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

/**
 * Calculates complex roots of unity.
 *
 * @param N The number of roots to calculate
 * @return A list of all complex numbers \f$z\f$ such that \f$z^N = 1\f$.
 */
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


/**
 * Load a set of filter coefficients
 *
 * @param filtername String specifying a filter. There should be a corresponding
 *                file in the RUNTIME_DIR called <filtername>.dat.
 * @param type    Either ANALYSIS_FILTER or SYNTHESIS_FILTER
 * @param nchans  The number of channels that this filter will be applied to.
 *                For both ANALYSIS and SYNTHESIS filters, this should be
 *                the number of ANALYSIS channels.
 *
 * @todo Change nchans so that it uses `vm&rarr;nfine_chans`
 */
void vmLoadFilter( vcsbeam_context *vm, char *filtername, filter_type type, int nchans )
{
    // Make sure function arguments are not NULL
    if (filtername == NULL)
    {
        fprintf( stderr, "error: vmLoadFilter: "
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
        fprintf( stderr, "error: vmLoadFilter: "
                "could not find filter '%s'\n", filtername );
        exit(EXIT_FAILURE);
    }

    // Allocate new memory and initialise variables
    // If the specified kind of filter currently already exists, override it
    pfb_filter *filter;
    if (type == ANALYSIS_FILTER)
    {
        if (vm->analysis_filter != NULL)
            free_pfb_filter( vm->analysis_filter );

        vm->analysis_filter = (pfb_filter *)malloc( sizeof(pfb_filter) );
        filter = vm->analysis_filter;
    }
    else // if (type == SYNTHESIS_FILTER)
    {
        if (vm->synth_filter != NULL)
            free_pfb_filter( vm->synth_filter );

        vm->synth_filter = (pfb_filter *)malloc( sizeof(pfb_filter) );
        filter = vm->synth_filter;
    }

    // Initialise variables
    filter->nchans  = nchans;
    filter->ncoeffs = 0;
    filter->type    = type;

    // Go through file and count the number of coefficients
    // We assume that the file contains ONLY ASCII floating point
    // numbers separated by whitespace.
    double dummy;
    int num_read;

    while ((num_read = fscanf( f, "%lf", &dummy )) != EOF)
        filter->ncoeffs++;

    // Check to see that there is at least 1 coefficient!
    if (filter->ncoeffs == 0)
    {
        fprintf( stderr, "error: vmLoadFilter: "
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
        fprintf( stderr, "warning: vmLoadFilter: "
                "number of channels (%d) does not divide evenly into the "
                "number of coefficients (%d)\n", nchans, filter->ncoeffs );
    }
    filter->ntaps = filter->ncoeffs / nchans;

    // Pre-calculate the twiddle factors
    filter->twiddles = roots_of_unity( nchans );

    // Close the file
    fclose( f );
}

/**
 * Scale the filter coefficients by constant factor.
 *
 * @param type         Either ANALYSIS_FILTER or SYNTHESIS_FILTER
 * @param scale_factor The scaling factor
 */
void vmScaleFilterCoeffs( vcsbeam_context *vm, filter_type type, double scale_factor )
{
    // Get the requested filter
    pfb_filter *filter;
    if (type == ANALYSIS_FILTER)
        filter = vm->analysis_filter;
    else // if (type == SYNTHESIS_FILTER)
        filter = vm->synth_filter;

    // If there is no such filter, silently do nothing
    if (filter == NULL)
        return;

    // Loop through and multiply the scale factor to each coefficient
    int i;
    for (i = 0; i < filter->ncoeffs; i++)
        filter->coeffs[i] *= scale_factor;
}

/**
 * Free the memory associated with a PFB filter.
 *
 * @param filter The filter to be freed/
 *
 * This function frees both the member variables and the filter itself.
 */
void free_pfb_filter( pfb_filter *filter )
/* Free the memory allocated in vmLoadFilter()
 */
{
    free( filter->coeffs );
    free( filter );

    filter = NULL;
}
