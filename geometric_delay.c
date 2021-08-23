/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <mwalib.h>
#include <star/pal.h>
#include <star/palmac.h>
#include "geometric_delay.h"
#include "beam_common.h"

void calc_geometric_delays(
        geometric_delays  *gdelays,
        struct beam_geom  *beam_geom_vals )
/* Calculate the geometric delay (in radians) for the given pointings
 */
{
    uintptr_t p, a, c; // Counters for (p)ointing, (a)ntenna, (c)hannel

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
    double L, w, Delta_t, phase;

    for (p = 0; p < gdelays->npointings; p++)
    {
        for (a = 0; a < gdelays->nant; a++)
        {
            // Get the location and cable length for this antenna
            E     = gdelays->obs_metadata->antennas[a].east_m;
            N     = gdelays->obs_metadata->antennas[a].north_m;
            H     = gdelays->obs_metadata->antennas[a].height_m;
            cable = gdelays->obs_metadata->antennas[a].electrical_length_m;
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

            for (c = 0; c < gdelays->nchan; c++)
            {
                // Eq. (3) in Ord et al. (2019)
                // NB: Again, the sign flip (compared to the paper) is
                // unexplained.
                phase = -2.0 * M_PI * Delta_t * gdelays->chan_freqs_hz[c];
                gdelays->phi[PHI_IDX(p, a, c, gdelays->nant, gdelays->nchan)] =
                    make_cuDoubleComplex( cos( phase ), sin( phase ));
            }
        }
    }
}

void create_geometric_delays(
        geometric_delays  *gdelays,
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        uintptr_t          coarse_chan_idx,
        uintptr_t          npointings )
/* Allocates memory for the geometric delay arrays ("phi") on both host and device.
 * Free with free_geometric_delays()
 */
{
    gdelays->npointings   = npointings;
    gdelays->nant         = obs_metadata->num_ants;
    gdelays->nchan        = vcs_metadata->num_fine_chans_per_coarse;
    gdelays->obs_metadata = obs_metadata;

    // Get a pointer to the array of fine channel frequencies
    // This is a (contiguous) subset of an array already provided
    // by mwalib, but we need to jump into it at the right coarse channel
    gdelays->chan_freqs_hz = &(obs_metadata->metafits_fine_chan_freqs[coarse_chan_idx*gdelays->nchan]);

    // Allocate memory
    size_t size = gdelays->npointings * gdelays->nant * gdelays->nchan * sizeof(cuDoubleComplex);

    cudaMallocHost( (void **)&(gdelays->phi), size );
    cudaCheckErrors( "error: create_geometric_delays: cudaMallocHost failed" );

    cudaMalloc( (void **)&(gdelays->d_phi), size );
    cudaCheckErrors( "error: create_geometric_delays: cudaMalloc failed" );
}

void free_geometric_delays( geometric_delays *gdelays )
/* Free memory allocated with create_geometric_delays()
 */
{
    cudaFreeHost( gdelays->phi );
    cudaCheckErrors( "(free_geometric_delays) cudaFreeHost failed" );

    cudaFree( gdelays->d_phi );
    cudaCheckErrors( "(free_geometric_delays) cudaFree failed" );
}

void push_geometric_delays_to_device( geometric_delays *gdelays )
/* Copy host memory block to device
 */
{
    size_t size = gdelays->npointings * gdelays->nant * gdelays->nchan * sizeof(cuDoubleComplex);
    cudaMemcpyAsync( gdelays->d_phi, gdelays->phi, size, cudaMemcpyHostToDevice, 0 );
    cudaCheckErrors( "error: push_geometric_delays_to_device: cudaMemcpyAsync failed" );
}


void calc_beam_geom(
        double            ras_hours,
        double            decs_degs,
        double            mjd,
        struct beam_geom  *bg )
{
    // Calculate geometry of pointings

    double intmjd;
    double fracmjd;
    double lmst;
    double mean_ra, mean_dec, ha;
    double az, el;

    double unit_N;
    double unit_E;
    double unit_H;

    double pr = 0, pd = 0, px = 0, rv = 0, eq = 2000, ra_ap = 0, dec_ap = 0;

    /* get mjd */
    intmjd = floor(mjd);
    fracmjd = mjd - intmjd;

    /* get requested Az/El from command line */
    mjd2lst( mjd, &lmst );

    /* for the look direction <not the tile> */

    mean_ra = ras_hours * PAL__DH2R;
    mean_dec = decs_degs * PAL__DD2R;

    palMap(mean_ra, mean_dec, pr, pd, px, rv, eq, mjd, &ra_ap, &dec_ap);

    // Lets go mean to apparent precess from J2000.0 to EPOCH of date.

    ha = palRanorm( lmst - ra_ap ); // in radians

    /* now HA/Dec to Az/El */

    palDe2h( ha, dec_ap, MWA_LATITUDE_RADIANS, &az, &el );

    /* now we need the direction cosines */

    unit_N = cos(el) * cos(az);
    unit_E = cos(el) * sin(az);
    unit_H = sin(el);

    // Populate a structure with some of the calculated values
    bg->mean_ra  = mean_ra;
    bg->mean_dec = mean_dec;
    bg->az       = az;
    bg->el       = el;
    bg->lmst     = lmst;
    bg->fracmjd  = fracmjd;
    bg->intmjd   = intmjd;
    bg->unit_N   = unit_N;
    bg->unit_E   = unit_E;
    bg->unit_H   = unit_H;
}

