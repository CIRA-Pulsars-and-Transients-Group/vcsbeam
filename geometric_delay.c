/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mwalib.h>
#include "geometric_delay.h"
#include "beam_common.h"

void calc_geometric_delays(
        double ***phi,
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        uintptr_t          coarse_chan_idx,
        struct beam_geom  *beam_geom_vals,
        uintptr_t          npointings )
/* Calculate the geometric delay (in radians) for the given
 *   pointings, antennas, channels
 */
{
    uintptr_t nant = obs_metadata->num_ants;
    uintptr_t nchan = vcs_metadata->num_fine_chans_per_coarse;
    uintptr_t p, a, c; // Counters for (p)ointing, (a)ntenna, (c)hannel

    // Get a pointer to the array of fine channel frequencies
    // This is a (contiguous) subset of an array already provided
    // by mwalib, but we need to jump into it at the right coarse channel
    uintptr_t coarse_chan = vcs_metadata->provided_coarse_chan_indices[coarse_chan_idx];
    double *fine_chan_freqs_hz = &(obs_metadata->metafits_fine_chan_freqs[coarse_chan*nchan]);

    double E, N, H; // Location of the antenna: (E)ast, (N)orth, (H)eight
    double cable;   // The cable length for a given antenna

    // Reference location (point on the ground to phase up to) and reference
    // cable length.
    // These are somewhat arbitrary -- the only point is to choose values such
    // that the possibility that a whole number of sample delays are required
    // is minimised. The "centre" of the array seems like a
    // good candidate.
    double Eref      = 0.0;
    double Nref      = 0.0;
    double Href      = MWA_ALTITUDE_METRES;
    double cable_ref = 0.0;

    // Other various intermediate products
    double L, w, Delta_t;

    for (p = 0; p < npointings; p++)
    {
        for (a = 0; a < nant; a++)
        {
            // Get the location and cable length for this antenna
            E     = obs_metadata->antennas[a].east_m;
            N     = obs_metadata->antennas[a].north_m;
            H     = obs_metadata->antennas[a].height_m;
            cable = obs_metadata->antennas[a].electrical_length_m;
            L     = cable - cable_ref;

            // Eq. (1) in Ord et al. (2019)
            // Get the geometric delay associated with the given pointing.
            // This equation is actually equivalent to Eq. (1) in Ord et al.
            // (2019), even though it is expressed somewhat differently (i.e.
            // as a dot product with the pre-calculated direction vector).
            w = (E - Eref)*beam_geom_vals[p].unit_E +
                (N - Nref)*beam_geom_vals[p].unit_N +
                (H - Href)*beam_geom_vals[p].unit_H;

            // Eq. (2) in Ord et al. (2019)
            // NB: The sign is different compared to the paper. I haven't
            // decided if it's just down to conventions, or whether it's a
            // mistake in the paper. In any case, a minus here gives the
            // correct answer.
            Delta_t = (w - L)/SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;

            for (c = 0; c < nchan; c++)
            {
                // Eq. (3) in Ord et al. (2019)
                // NB: Again, the sign flip (compared to the paper) is
                // unexplained.
                phi[p][a][c] = -2.0 * M_PI * Delta_t * fine_chan_freqs_hz[c];
            }
        }
    }
}

double ***malloc_geometric_delays(
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        uintptr_t          npointings )
/* Allocates memory for the geometric delay arrays ("phi").
 * Free with free_geometric_delays().
 */
{
    uintptr_t nant = obs_metadata->num_ants;
    uintptr_t nchan = vcs_metadata->num_fine_chans_per_coarse;
    uintptr_t p, a; // Counters for (p)ointing, (a)ntenna

    double ***phi = (double ***)malloc( npointings * sizeof(double **) );
    for (p = 0; p < npointings; p++)
    {
        phi[p] = (double **)malloc( nant * sizeof(double *) );
        for (a = 0; a < nant; a++)
        {
            phi[p][a] = (double *)malloc( nchan * sizeof(double) );
        }
    }

    return phi;
}


void free_geometric_delays(
        double           ***phi,
        MetafitsMetadata   *obs_metadata,
        uintptr_t           npointings )
/* Free memory created with malloc_geometric_delays()
 */
{
    uintptr_t nant = obs_metadata->num_ants;
    uintptr_t p, a; // Counters for (p)ointing, (a)ntenna

    for (p = 0; p < npointings; p++)
    {
        for (a = 0; a < nant; a++)
            free( phi[p][a] );
        free( phi[p] );
    }
    free( phi );
}

