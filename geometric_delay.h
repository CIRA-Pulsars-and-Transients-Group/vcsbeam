/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __GEOMETRIC_DELAY_H__
#define __GEOMETRIC_DELAY_H__

#include <cuComplex.h>
#include <mwalib.h>
#include "beam_common.h"

/* In order to communicate easily with the GPU, the "phi" array, which
 * contains the complex-valued geometric delay terms, is implemented
 * as a 1D array, which uses the following indexing macro to access
 * the correct term for a given (p)ointing, (a)ntenna, and (c)hannel.
 */

#define PHI_IDX(p,a,c,na,nc)   ((p) * (nc)*(na)  + \
                                (a) * (nc)       + \
                                (c))

typedef struct geometric_delays_t {
    cuDoubleComplex   *phi;
    cuDoubleComplex   *d_phi;
    uintptr_t          npointings;
    uintptr_t          nant;
    uintptr_t          nchan;
    double            *chan_freqs_hz;
    MetafitsMetadata  *obs_metadata;
} geometric_delays;

/* Calculate the geometric delay (in radians) for the given pointings
 */
void calc_geometric_delays(
        geometric_delays  *gdelays,
        struct beam_geom  *beam_geom_vals );

/* Allocates memory for the geometric delay arrays ("phi") on both host and device.
 * Free with free_geometric_delays()
 */
void create_geometric_delays(
        geometric_delays  *gdelays,
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        uintptr_t          coarse_chan,
        uintptr_t          npointings );


/* Free memory allocated with create_geometric_delays()
 */
void free_geometric_delays( geometric_delays *gdelays );

/* Copy host memory block to device
 */
void push_geometric_delays_to_device( geometric_delays *gdelays );

#endif
