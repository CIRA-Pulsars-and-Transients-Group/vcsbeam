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
#include "geometry.h"
#include "jones.h"

void calc_geometric_delays(
        struct beam_geom  *beam_geom_vals,
        uint32_t           freq_hz,
        MetafitsMetadata  *obs_metadata,
        cuDoubleComplex   *phi )
{
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

    for (uintptr_t a = 0; a < obs_metadata->num_ants; a++)
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
        w = (E - Eref)*beam_geom_vals->unit_E +
            (N - Nref)*beam_geom_vals->unit_N +
            (H - Href)*beam_geom_vals->unit_H;

        // Eq. (2) in Ord et al. (2019)
        // NB: The sign is different compared to the paper. I haven't
        // decided if it's just down to conventions, or whether it's a
        // mistake in the paper. In any case, a minus here gives the
        // correct answer.
        Delta_t = (w - L)/SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;

        // Eq. (3) in Ord et al. (2019)
        // NB: Again, the sign flip (compared to the paper) is
        // unexplained.
        phase = -2.0 * M_PI * Delta_t * freq_hz;
        phi[a] = make_cuDoubleComplex( cos( phase ), sin( phase ));
    }
}

void calc_all_geometric_delays(
        geometric_delays  *gdelays,
        struct beam_geom  *beam_geom_vals )
/* Calculate the geometric delay (in radians) for the given pointings
 */
{
    cuDoubleComplex phi[gdelays->nant]; // <-- Temporary location for result

    for (uintptr_t p = 0; p < gdelays->npointings; p++)
    {
        for (uintptr_t c = 0; c < gdelays->nchan; c++)
        {
            calc_geometric_delays(
                    &beam_geom_vals[p],
                    gdelays->chan_freqs_hz[c],
                    gdelays->obs_metadata,
                    phi );

            for (uintptr_t a = 0; a < gdelays->nant; a++)
                gdelays->phi[PHI_IDX(p, a, c, gdelays->nant, gdelays->nchan)] = phi[a];
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



void mjd2lst(double mjd, double *lst)
{
    // Greenwich Mean Sidereal Time to LMST
    // east longitude in hours at the epoch of the MJD
    double lmst = palRanorm(palGmst(mjd) + MWA_LONGITUDE_RADIANS);

    *lst = lmst;
}


double parse_ra( char* ra_hhmmss )
/* Parse a string containing a right ascension in hh:mm:ss format into
 * a double in units of hours
 */
{
    if (ra_hhmmss == NULL)
        return 0.0;

    int ih=0, im=0, J=0;
    double fs=0., ra_rad=0.;

    sscanf(ra_hhmmss, "%d:%d:%lf", &ih, &im, &fs);

    palDtf2r(ih, im, fs, &ra_rad, &J);

    if (J != 0) { // pal returned an error
        fprintf(stderr,"Error parsing %s as hhmmss\npal error code: j=%d\n",ra_hhmmss,J);
        fprintf(stderr,"ih = %d, im = %d, fs = %lf\n", ih, im, fs);
        exit(EXIT_FAILURE);
    }

    return ra_rad*PAL__DR2H;
}

double parse_dec( char* dec_ddmmss )
/* Parse a string containing a declination in dd:mm:ss format into
 * a double in units of degrees
 */
{
    if (dec_ddmmss == NULL)
        return 0.0;

    int id=0, im=0, J=0, sign=0;
    double fs=0., dec_rad=0.;
    char id_str[16];

    sscanf(dec_ddmmss, "%s:%d:%lf", id_str, &im, &fs);

    if (id_str[0] == '-') {
        sign = -1;
    }
    else {
        sign = 1;
    }
    sscanf(dec_ddmmss, "%d:%d:%lf", &id, &im, &fs);
    id = id*sign;
    palDaf2r(id, im, fs, &dec_rad, &J);

    if (J != 0) {
        fprintf(stderr,"Error parsing %s as dd:mm:ss - got %d:%d:%f -- error code %d\n",dec_ddmmss,id,im,fs,J);
        exit(EXIT_FAILURE);
    }

    return dec_rad*PAL__DR2D*sign;
}


double calc_array_factor(
        MetafitsMetadata *obs_metadata,
        uint32_t          freq_hz,
        struct beam_geom *bg1,
        struct beam_geom *bg2 )
/* An implementation of Eq. (A4) in Meyers et al. (2017)
 */
{
    double array_factor;

    uintptr_t nant = obs_metadata->num_ants;
    cuDoubleComplex w[nant], psi[nant];

    calc_geometric_delays( bg1, freq_hz, obs_metadata, w );
    calc_geometric_delays( bg2, freq_hz, obs_metadata, psi );

    cuDoubleComplex cumsum = make_cuDoubleComplex( 0.0, 0.0 );
    for (uintptr_t a = 0; a < nant; a++)
    {
        cumsum = cuCadd( cumsum, cuCmul( cuConj(w[a]), psi[a] ) );
    }

    array_factor = cuCabs( cumsum );
    array_factor *= array_factor / nant;

    return array_factor;
}
