/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __GEOMETRIC_DELAY_H__
#define __GEOMETRIC_DELAY_H__

#include <mwalib.h>

/* Calculate the geometric delay (in radians) for the given
 *   pointings, antennas, channels
 */
void calc_geometric_delays(
        double ***phi,
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        uintptr_t          coarse_chan_idx,
        struct beam_geom  *beam_geom_vals,
        uintptr_t          npointings );

/* Allocates memory for the geometric delay arrays ("phi").
 * Free with free_geometric_delays().
 */
double ***malloc_geometric_delays(
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        uintptr_t          npointings );


/* Free memory created with malloc_geometric_delays()
 */
void free_geometric_delays(
        double           ***phi,
        MetafitsMetadata   *obs_metadata,
        uintptr_t           npointings );

#endif
