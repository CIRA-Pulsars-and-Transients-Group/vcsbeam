/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <cuComplex.h>
#include "star/pal.h"
#include "star/palmac.h"
#include "psrfits.h"
#include "fitsio.h"
#include <string.h>
#include <mwalib.h>
#include "beam_common.h"
#include "mwa_hyperbeam.h"
#include "calibration.h"

/* make a connection to the MWA database and get the antenna positions.
 * Then: calculate the geometric delay to a source for each antenna
 *
 * Initially geocentric to match as best as possible the output of calc.
 *
 */


/* This now also reads in RTS calibration files and generates calibration matrices
 */

#define MAXREQUEST 3000000

//=====================//

void create_delays_amps_from_metafits( MetafitsMetadata *obs_metadata, uint32_t ***delays, double ***amps )
/* Costruct "delay" and "amp" arrays required by Hyperbeam.
 * Delays can be brought directly across from the MetafitsMetadata struct.
 * The "conversion" from delays to amps is straightforward: All amps are "1.0"
 * unless the delay is 32, in which case the amp should be "0.0".
 * (The value 32 is the number assigned to a faulty or broken dipole.)
 *
 * Both the delays array (input) and the amps array (output) should
 * have dimensions [ninputs][ndipoles]
 *
 * This function allocates memory for both delays and amps. This can be freed
 * by a call to free_delays_amps().
 */
{
    int ndipoles = obs_metadata->num_delays;
    int ninputs  = obs_metadata->num_rf_inputs;

    *delays = (uint32_t **)malloc( ninputs * sizeof(uint32_t *) );
    *amps = (double **)malloc( ninputs * sizeof(double *) );

    int input, dipole;

    // For the delays, we will just point each delay array to the
    // corresponding array in the metafits, but for amps, we need
    // to allocate new memory for each set of ndipoles (=16) numbers
    // and populate them
    for (input = 0; input < ninputs; input++)
    {
        (*delays)[input] = obs_metadata->rf_inputs[input].dipole_delays;

        (*amps)[input] = (double *)malloc( ndipoles * sizeof(double) );
        for (dipole = 0; dipole < ndipoles; dipole++)
        {
            (*amps)[input][dipole] = ((*delays)[input][dipole] == 32 ? 0.0 : 1.0);
        }
    }
}

void free_delays_amps( MetafitsMetadata *obs_metadata, uint32_t **delays, double **amps )
/* Frees the memory allocated in create_delays_amps_from_metafits().
 */
{
    int ninputs  = obs_metadata->num_rf_inputs;
    int input;
    for (input = 0; input < ninputs; input++)
        free( amps[input] );

    free( amps );
    free( delays );
}

void create_antenna_lists( MetafitsMetadata *obs_metadata, uint32_t *polX_idxs, uint32_t *polY_idxs )
/* Creates a list of indexes into the data for the X and Y polarisations,
 * ordered by antenna number. Assumption: polX_idxs and polY_idxs have
 * sufficient allocated memory.
 */
{
    // Go through the rf_inputs and construct the lookup table for the antennas
    unsigned int ninputs = obs_metadata->num_rf_inputs;
    unsigned int i, ant;
    for (i = 0; i < ninputs; i++)
    {
        ant = obs_metadata->rf_inputs[i].ant;
        if (*(obs_metadata->rf_inputs[i].pol) == 'X')
            polX_idxs[ant] = i;
        else // if (*(obs_metadata->rf_inputs.pol) == 'Y')
            polY_idxs[ant] = i;
    }
}

int hash_dipole_config( double *amps )
/* In order to avoid recalculating the FEE beam for repeated dipole
 * configurations, we have to keep track of which configurations have already
 * been calculated. We do this through a boolean array, and this function
 * converts dipole configurations into indices of this array. In other words,
 * this function _assigns meaning_ to the array.
 *
 * Since dead dipoles are relatively rare, we only consider configurations
 * in which up to two dipoles are dead. Any more than that and the we can
 * recalculate the Jones matrix with minimal entropy. In this case, this
 * function returns -1. The other remaining cases are:
 *
 *   0  dead dipoles = 1   configuration
 *   1  dead dipole  = 16  configurations
 *   2  dead dipoles = 120 configurations
 *   16 dead dipoles = 1   configuration
 *
 * for a grand total of 138 indices. They are ordered as follows:
 *
 *  idx  configuration (0=dead, 1=live)
 *   0   [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   1   [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   2   [1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   3   [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   ...
 *   16  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
 *   17  [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   18  [0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   19  [0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   ...
 *   31  [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
 *   32  [1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   32  [1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   ...
 *   136 [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0]
 *   136 [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
 *   137 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 *
 *   This function defines "dead" as amps=0.0, "live" otherwise.
 */
{
    // The return value = the "hashed" array index
    int idx;
    int d1 = 0, d2 = 0; // The "locations" of the (up to) two dead dipoles

    // Count the number of dead dipoles
    int ndead = 0;
    int ndipoles = 16;
    int i;
    for (i = 0; i < ndipoles; i++)
        ndead += (amps[i] == 0.0);

    // Convert configuration into an array index (the "hash" function)
    switch (ndead)
    {
        case 0:
            idx = 0;
            break;
        case 1:
            for (i = 0; i < ndipoles; i++)
            {
                if (amps[i] == 0.0)
                {
                    d1 = i;
                    break;
                }
            }
            idx = d1 + 1; // The "1-dead-dipole" configs start at idx 1
            break;
        case 2:
            // Find the locations of the two dead dipoles
            d1 = -1;
            for (i = 0; i < ndipoles; i++)
            {
                if (amps[i] == 0.0)
                {
                    if (d1 == -1)
                        d1 = i;
                    else
                    {
                        d2 = i;
                        break;
                    }
                }
            }
            // The hash function for two dead dipoles
            // (The "2-dead-dipole" configs start at idx 17
            idx = 16*d1 - ((d1 + 2)*(d1 + 1))/2 + d2 + 17;
            break;
        case 16:
            idx = 137;
            break;
        default: // any configuration with >2 dead dipoles
            idx = -1;
            break;
    }

    return idx;
}

double parse_dec( char* dec_ddmmss ) {
/* Parse a string containing a declination in dd:mm:ss format into
 * a double in units of degrees
 */

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

double parse_ra( char* ra_hhmmss ) {
/* Parse a string containing a right ascension in hh:mm:ss format into
 * a double in units of hours
 */

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


void mjd2lst(double mjd, double *lst) {

    // Greenwich Mean Sidereal Time to LMST
    // east longitude in hours at the epoch of the MJD
    double lmst = palRanorm(palGmst(mjd) + MWA_LONGITUDE_RADIANS);

    *lst = lmst;
}

void get_jones(
        // an array of pointings [pointing][ra/dec][characters]
        int                    npointing, // number of pointings
        VoltageMetadata       *vcs_metadata,
        MetafitsMetadata      *obs_metadata,
        int                    coarse_chan_idx,
        struct                 calibration *cal,
        cuDoubleComplex     ***D,
        FEEBeam               *beam,
        uint32_t             **delays,
        double               **amps,
        struct beam_geom       beam_geom_vals[],
        cuDoubleComplex    ****complex_weights_array,  // output: cmplx[npointing][ant][ch][pol]
        cuDoubleComplex    ****invJi )                 // output: invJi[ant][ch][pol][pol]
{

    // Give "shorthand" variables for often-used values in metafits
    int coarse_chan    = vcs_metadata->common_coarse_chan_indices[coarse_chan_idx];
    long int frequency = obs_metadata->metafits_coarse_chans[coarse_chan].chan_start_hz;
    //int nant           = obs_metadata->num_ants;
    int nchan          = obs_metadata->num_volt_fine_chans_per_coarse;
    int npol           = obs_metadata->num_ant_pols;   // (X,Y)
    int chan_width     = obs_metadata->corr_fine_chan_width_hz;
    int ninput         = obs_metadata->num_rf_inputs;
    unsigned int samples_per_sec = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;

    int rf_input;     // For counting through nstation*npol rows in the metafits file
    int ant;     // Antenna number
    int pol;     // Polarisation number
    int ch;      // Channel number
    int p1, p2;  // Counters for polarisation

    /* easy -- now the positions from the database */

    double phase;

    cuDoubleComplex E[npol*npol];               // Model Jones in Desired Direction
    cuDoubleComplex Ji[npol*npol];              // Gain in Desired Direction
    double        P[npol*npol];               // Parallactic angle correction rotation matrix

    // Choose a reference tile
    int refinp = 84; // Tile012 (?)
    Antenna ref_ant = obs_metadata->antennas[refinp];
    double N_ref = ref_ant.north_m;
    double E_ref = ref_ant.east_m;
    double H_ref = ref_ant.height_m;

    int n;

    long int freq_ch;

    double cable;

    double integer_phase;
    double geometry, delay_time, delay_samples, cycles_per_sample;

    int nconfigs = 139;
    const int multiple_dead = nconfigs - 1; // This indexn is reserved for configurations
                                            // with three or more dead dipoles
    int config_idx;
    double *jones[nconfigs]; // (see hash_dipole_configs() for explanation of this array)

    double uv_angle;
    cuDoubleComplex uv_phase; // For the UV phase correction

    double Fnorm;

    // Set settings for the FEE2016 beam model using Hyperbeam
    int zenith_norm = 1; // Boolean value: unsure if/how this makes a difference

    for (int p = 0; p < npointing; p++)
    {
        // Reset the Jones matrices (for the FEE beam)
        for (n = 0; n < nconfigs; n++)
            jones[n] = NULL; // i.e. no Jones matrices have been calculated for any configurations so far

        parallactic_angle_correction_fee2016(
                P,                       // output = rotation matrix
                MWA_LATITUDE_RADIANS,    // observing latitude (radians)
                beam_geom_vals[p].az, (PAL__DPIBY2-beam_geom_vals[p].el));   // azimuth & zenith angle of pencil beam

        // Everything from this point on is frequency-dependent
        for (ch = 0; ch < nchan; ch++) {

            // Calculating direction-dependent matrices
            freq_ch = frequency + ch*chan_width;    // The frequency of this fine channel

            // Calculate the UV phase correction for this channel
            uv_angle = cal->phase_slope*freq_ch + cal->phase_offset;
            uv_phase = make_cuDoubleComplex( cos(uv_angle), sin(uv_angle) );

            for (rf_input = 0; rf_input < (int)(ninput); rf_input++) {

                // Get the antenna and polarisation number from the rf_input
                ant         = obs_metadata->rf_inputs[rf_input].ant;
                pol         = *(obs_metadata->rf_inputs[rf_input].pol) - 'X'; // 'X' --> 0; 'Y' --> 1

                // FEE2016 beam:
                // Check to see whether or not this configuration has already been calculated.
                // The point of this is to save recalculating the jones matrix, which is
                // computationally expensive.
                config_idx = hash_dipole_config( amps[rf_input] );
                if (config_idx == -1)
                    config_idx = multiple_dead;

                if (ch == 0 && (jones[config_idx] == NULL || config_idx == multiple_dead))
                {
                    // The Jones matrix for this configuration has not yet been calculated, so do it now.
                    // The FEE beam only needs to be calculated once per coarse channel, because it will
                    // not produce unique answers for different fine channels within a coarse channel anyway
                    // (it only calculates the jones matrix for the nearest coarse channel centre)
                    // Strictly speaking, the condition (ch == 0) above is redundant, as the dipole configuration
                    // array takes care of that implicitly, but I'll leave it here so that the above argument
                    // is "explicit" in the code.
                    jones[config_idx] = calc_jones( beam, beam_geom_vals[p].az, PAL__DPIBY2 - beam_geom_vals[p].el, frequency + chan_width/2,
                            (unsigned int*)delays[rf_input], amps[rf_input], zenith_norm );
                }

                // "Convert" the real jones[8] output array into out complex E[4] matrix
                for (n = 0; n < npol*npol; n++)
                {
                    E[n] = make_cuDoubleComplex(jones[config_idx][n*2], jones[config_idx][n*2+1]);
                }

                // Apply parallactic angle correction if Hyperbeam was used
                mult2x2d_RxC( P, E, E );  // Ji = P x Ji (where 'x' is matrix multiplication)

                mult2x2d(D[ant][ch], E, Ji); // the gain in the desired look direction

                // Apply the UV phase correction (to the bottom row of the Jones matrix)
                Ji[2] = cuCmul( Ji[2], uv_phase );
                Ji[3] = cuCmul( Ji[3], uv_phase );

                // Calculate the complex weights array
                if (complex_weights_array != NULL)
                {
                    cable = obs_metadata->rf_inputs[rf_input].electrical_length_m - ref_ant.electrical_length_m;
                    double El = obs_metadata->rf_inputs[rf_input].east_m;
                    double N  = obs_metadata->rf_inputs[rf_input].north_m;
                    double H  = obs_metadata->rf_inputs[rf_input].height_m;

                    // shift the origin of ENH to Antenna 0 and hoping the Far Field Assumption still applies ...

                    geometry = (El - E_ref)*beam_geom_vals[p].unit_E +
                               (N  - N_ref)*beam_geom_vals[p].unit_N +
                               (H  - H_ref)*beam_geom_vals[p].unit_H;

                    delay_time = (geometry - cable)/(SPEED_OF_LIGHT_IN_VACUUM_M_PER_S);
                    delay_samples = delay_time * samples_per_sec;

                    // freq should be in cycles per sample and delay in samples
                    // which essentially means that the samples_per_sec cancels

                    // and we only need the fractional part of the turn
                    cycles_per_sample = (double)freq_ch/samples_per_sec;

                    phase = cycles_per_sample*delay_samples;
                    phase = modf( phase, &integer_phase );

                    phase *= -2.0*M_PI;

                    // Store result for later use
                    complex_weights_array[p][ant][ch][pol] =
                        make_cuDoubleComplex( cos( phase ), sin( phase ));
                }

                // Now, calculate the inverse Jones matrix
                if (invJi != NULL)
                {
                    if (pol == 0) { // This is just to avoid doing the same calculation twice

                        // Ord's original comment for the following line is:
                        // "The RTS conjugates the sky so beware"
                        conj2x2( Ji, Ji );

                        Fnorm = norm2x2( Ji, Ji );

                        if (Fnorm != 0.0)
                            inv2x2S( Ji, invJi[ant][ch] );
                        else {
                            for (p1 = 0; p1 < npol;  p1++)
                            for (p2 = 0; p2 < npol;  p2++)
                                invJi[ant][ch][p1][p2] = make_cuDoubleComplex( 0.0, 0.0 );
                        }
                    }
                }

            } // end loop through antenna/pol (rf_input)
        } // end loop through fine channels (ch)

        // Free Jones matrices from hyperbeam -- in prep for reclaculating the next pointing
        for (n = 0; n < nconfigs; n++)
        {
            if (jones[n] != NULL)
                free( jones[n] );
        }
    } // end loop through pointings (p)

}

void parallactic_angle_correction_fee2016(
    double *P,    // output rotation matrix
    double lat,   // observing latitude (radians)
    double az,    // azimuth angle (radians)
    double za)    // zenith angle (radians)
{
    // The FEE beam Jones matrix is
    //   [ Qθ  Pθ ]
    //   [ Qφ  Pφ ]
    // RTS calibration solution is
    //   [ PP  PQ ]
    //   [ QP  QQ ]
    // Therefore, the parallactic rotation must be
    //   [  0   1 ] [ cos(χ)  -sin(χ) ]  =  [  sin(χ)  cos(χ) ]
    //   [  1   0 ] [ sin(χ)   cos(χ) ]     [  cos(χ) -sin(χ) ]

    double el = PAL__DPIBY2 - za;

    double ha, dec;
    palDh2e(az, el, lat, &ha, &dec);
    double pa = palPa( ha, dec, lat );

    P[0] = sin(pa);
    P[1] = cos(pa);
    P[2] = cos(pa);
    P[3] = -sin(pa);
}

void calc_beam_geom(
        char              pointing_array[][2][64],
        int               npointing,
        double            mjd,
        struct beam_geom  bg[] )
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

    double dec_degs;
    double ra_hours;

    double pr = 0, pd = 0, px = 0, rv = 0, eq = 2000, ra_ap = 0, dec_ap = 0;

    /* get mjd */
    intmjd = floor(mjd);
    fracmjd = mjd - intmjd;

    /* get requested Az/El from command line */
    mjd2lst( mjd, &lmst );

    for (int p = 0; p < npointing; p++)
    {
        dec_degs = parse_dec( pointing_array[p][1] );
        ra_hours = parse_ra( pointing_array[p][0] );

        /* for the look direction <not the tile> */

        mean_ra = ra_hours * PAL__DH2R;
        mean_dec = dec_degs * PAL__DD2R;

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
        bg[p].mean_ra  = mean_ra;
        bg[p].mean_dec = mean_dec;
        bg[p].az       = az;
        bg[p].el       = el;
        bg[p].lmst     = lmst;
        bg[p].fracmjd  = fracmjd;
        bg[p].intmjd   = intmjd;
        bg[p].unit_N   = unit_N;
        bg[p].unit_E   = unit_E;
        bg[p].unit_H   = unit_H;
    }
}
