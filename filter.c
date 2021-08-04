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
#include "filter.h"


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


