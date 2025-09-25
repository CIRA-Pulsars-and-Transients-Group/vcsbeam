/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
// #include <cuComplex.h>
// #include <cuda_runtime.h>

#include <mwalib.h>
#include <star/pal.h>
#include <star/palmac.h>

#include "vcsbeam.h"
#include "gpu_macros.h"

/**
 * Calculates the phase delays required for phasing up each antenna.
 *
 * @param[in]  beam_geom_vals Struct containing pointing information for the
 *                            requested beams.
 * @param[in]  freq_hz        The frequency in Hz
 * @param[in]  obs_metadata   The observation metadata
 * @param[out] phi            The calculated phase delays
 *
 * This function calculates the phase delays for the given look-direction
 * and frequency, for each antenna. This consists of both a "geometric delay"
 * component, related to the different times a planar wavefront coming from
 * a particular direction arrives at each antenna, and a "cable delay"
 * component, related to the physical lengths of the cables connecting the
 * antennas to the rest of the system.
 *
 * The equations implemented here are the first three equations in
 * [Ord et al. (2019)](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/abs/mwa-tiedarray-processing-i-calibration-and-beamformation/E9A7A9981AE9A935C9E08500CA6A1C1E).
 */
void calc_geometric_delays(
        beam_geom         *beam_geom_vals,
        uint32_t           freq_hz,
        MetafitsMetadata  *obs_metadata,
        gpuDoubleComplex   *phi )
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

    uintptr_t i;
    Rfinput *Rf;
    for (i = 0; i < obs_metadata->num_rf_inputs; i++)
    {
        // Get a shortcut to this Rf input
        Rf = &(obs_metadata->rf_inputs[i]);

        // Only need to do this once per antenna, so just do it for X pols
        if (*(Rf->pol) == 'Y')
            continue;

        // Get the location and cable length for this antenna
        E     = Rf->east_m;
        N     = Rf->north_m;
        H     = Rf->height_m;
        cable = Rf->electrical_length_m;
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
        phase = 2.0 * M_PI * Delta_t * freq_hz;
        phi[Rf->ant] = make_gpuDoubleComplex( cos( phase ), sin( phase ));
    }
}

/**
 * Calculates the phase delays for all relevant frequencies and pointings.
 *
 * Calls calc_geometric_delays() for all frequencies in the relevant coarse
 * channel, and for all requested pointings.
 *
 * @todo Incorporate the `beam_geom` struct into the `vm` struct.
 */
void vmCalcPhi(
        vcsbeam_context   *vm,
        beam_geom         *beam_geom_vals )
/* Calculate the geometric delay (in radians) for the given pointings
 */
{
    geometric_delays *gdelays = &vm->gdelays;
    gpuDoubleComplex phi[gdelays->nant]; // <-- Temporary location for result

    uintptr_t p, c, a;
    for (p = 0; p < gdelays->npointings; p++)
    {
        for (c = 0; c < gdelays->nchan; c++)
        {
            calc_geometric_delays(
                    &beam_geom_vals[p],
                    gdelays->chan_freqs_hz[c],
                    gdelays->obs_metadata,
                    phi );

            for (a = 0; a < gdelays->nant; a++)
                gdelays->phi[PHI_IDX(p, a, c, gdelays->nant, gdelays->nchan)] = phi[a];
        }
    }
}

/**
 * Allocates memory for the delay phase arrays on both host and device.
 *
 * Free with free_geometric_delays()
 */
void vmCreateGeometricDelays( vcsbeam_context *vm )
{
    vm->gdelays.npointings   = vm->npointing;
    vm->gdelays.nant         = vm->obs_metadata->num_ants;
    vm->gdelays.nchan        = vm->nfine_chan;
    vm->gdelays.obs_metadata = vm->obs_metadata;

    // Create an array of fine channel frequencies.
    // This must be worked out ourselves (as opposed to grabbing it from mwalib)
    // in case we have done the fine channelisation ourselves (e.g. via a fine PFB)
    vm->gdelays.chan_freqs_hz = (double *)malloc( vm->nfine_chan * sizeof(double) );

    int coarse_chan_idx      = vm->coarse_chan_idxs_to_process[0];
    CoarseChannel C          = vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx];
    uint32_t fine_chan_width = vm->obs_metadata->coarse_chan_width_hz / vm->nfine_chan;
    uint32_t start_hz        = C.chan_start_hz;

    int c;
    for (c = 0; c < vm->nfine_chan; c++)
    {
        vm->gdelays.chan_freqs_hz[c] = start_hz + fine_chan_width*c;
    }

    // Allocate memory
    size_t size = vm->gdelays.npointings * vm->gdelays.nant * vm->gdelays.nchan * sizeof(gpuDoubleComplex);

    gpuMallocHost( (void **)&(vm->gdelays.phi), size );

    gpuMalloc( (void **)&(vm->gdelays.d_phi), size );
}

/**
 * Frees memory for the delay phase arrays on both host and device.
 *
 * @todo Convert free_geometric_delays() into a "vm" function.
 */
void free_geometric_delays( geometric_delays *gdelays )
/* Free memory allocated with create_geometric_delays()
 */
{
    free( gdelays->chan_freqs_hz );

    gpuHostFree( gdelays->phi );

    gpuFree( gdelays->d_phi );
}

/**
 * Copies the delay phase arrays from CPU memory to GPU memory.
 */
void vmPushPhi( vcsbeam_context *vm )
/* Copy host memory block to device
 */
{
    size_t size = vm->gdelays.npointings * vm->gdelays.nant * vm->gdelays.nchan * sizeof(gpuDoubleComplex);
    gpuMemcpyAsync( vm->gdelays.d_phi, vm->gdelays.phi, size, gpuMemcpyHostToDevice, 0 );
}

/**
 * Populates a `beam_geom` struct with pointing information derived from a
 * given set of RAs, Decs, and MJDs.
 *
 * @param[in]  ras_hours An array of RAs (in decimal hours)
 * @param[in]  decs_degs An array of Decs (in decimal degrees)
 * @param[in]  mjd       The Modified Julian Date
 * @param[out] bg        The struct containing various geometric quantities
 *
 * The quantities which are calculated are
 * | Quantity       | Description                    |
 * | -------------- | ------------------------------ |
 * | `bg->mean_ra`  | The mean RA of the pointing    |
 * | `bg->mean_dec` | The mean Dec of the pointing   |
 * | `bg->az`       | The azimuth of the pointing    |
 * | `bg->el`       | The elevation of the pointing  |
 * | `bg->lmst`     | The local mean sidereal time   |
 * | `bg->fracmjd`  | The fractional part of the MJD |
 * | `bg->intmjd`   | The integer part of the MJD    |
 * | `bg->unit_N`   | The normalised projection of the look-direction onto local North |
 * | `bg->unit_E`   | The normalised projection of the look-direction onto local East  |
 * | `bg->unit_H`   | The normalised projection of the look-direction onto local "Up"  |
 *
 * @todo Put the table describing the beam_geom struct where it belongs: with
 *       the documentation for the beam_geom struct!
 */
void calc_beam_geom(
        double            ras_hours,
        double            decs_degs,
        double            mjd,
        beam_geom        *bg )
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

    mean_ra = ras_hours * H2R;
    mean_dec = decs_degs * D2R;

    palMap(mean_ra, mean_dec, pr, pd, px, rv, eq, mjd, &ra_ap, &dec_ap);

    // Lets go mean to apparent precess from J2000.0 to EPOCH of date.

    ha = palRanorm( lmst - ra_ap ); // in radians

    /* now HA/Dec to Az/El */

    palDe2h( ha, dec_ap, MWA_LATITUDE_RADIANS, &az, &el );
    // ^-- Returns "geographic azimuth" and "elevation" (see documentation)

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

/**
 * Converts a decimal RA to the format "HH:MM:SS.S".
 *
 * @param[out] out   A buffer for the output string
 * @param[in]  in    The decimal RA
 * @param      sflag Whether to include a leading '+'
 *
 * The `out` buffer must already be allocated.
 *
 * @todo Figure out whether dec2hms() is *actually* just for RAs, or whether
 *       it is more general than that (and, if so, change the name to reflect
 *       this).
 */
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


/**
 * Convert MJD to LST.
 *
 * @param[in]  mjd  The Modified Julian Date
 * @param[out] lst  The Local Sidereal Time
 *
 * @todo Consider removing mjd2lst(), since it consists only of a single
 *       call to a `pal` function.
 */
void mjd2lst(double mjd, double *lst)
{
    // Greenwich Mean Sidereal Time to LMST
    // east longitude in hours at the epoch of the MJD
    double lmst = palRanorm(palGmst(mjd) + MWA_LONGITUDE_RADIANS);

    *lst = lmst;
}


/**
 * Parse an RA string.
 *
 * @param[in] ra_hhmmss A string representing an RA in the HH:MM:SS.S format
 * @return The RA in decimal hours.
 */
double parse_ra( char* ra_hhmmss )
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

    return ra_rad*R2H;
}

/**
 * Parse a Dec string.
 *
 * @param[in] dec_ddmmss A string representing a Dec in the DD:MM:SS.S format
 * @return The Dec in decimal degrees.
 */
double parse_dec( char* dec_ddmmss )
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

    return dec_rad*R2D*sign;
}


/**
 * Calculates the array factor.
 *
 * @param[in] obs_metadata  The observation metadata
 * @param[in] freq_hz       The frequency in Hz
 * @param[in] bg1           The look-direction
 * @param[in] bg2           The direction in which to calculate the array
 *                          factor
 * @return The array factor
 *
 * This function computes the array factor as described in the appendix of
 * [Meyers et al. (2017)](https://iopscience.iop.org/article/10.3847/1538-4357/aa8bba).
 */
double calc_array_factor(
        MetafitsMetadata *obs_metadata,
        uint32_t          freq_hz,
        beam_geom        *bg1,
        beam_geom        *bg2 )
/* An implementation of Eq. (A4) in Meyers et al. (2017)
 */
{
    double array_factor;

    uintptr_t nant = obs_metadata->num_ants;
    gpuDoubleComplex w[nant], psi[nant];

    calc_geometric_delays( bg1, freq_hz, obs_metadata, w );
    calc_geometric_delays( bg2, freq_hz, obs_metadata, psi );

    gpuDoubleComplex cumsum = make_gpuDoubleComplex( 0.0, 0.0 );
    uintptr_t a;
    for (a = 0; a < nant; a++)
    {
        cumsum = gpuCadd( cumsum, gpuCmul( gpuConj(w[a]), psi[a] ) );
    }

    array_factor = gpuCabs( cumsum ) / nant;
    array_factor *= array_factor;

    return array_factor;
}
