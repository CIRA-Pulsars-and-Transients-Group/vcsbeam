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
#include <star/pal.h>
#include <star/palmac.h>
#include "psrfits.h"
#include "fitsio.h"
#include <string.h>
#include <mwalib.h>
#include "beam_common.h"
#include <mwa_hyperbeam.h>
#include "calibration.h"
#include "primary_beam.h"

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
        cuDoubleComplex     ***B,
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

    cuDoubleComplex Ji[npol*npol];              // Gain in Desired Direction

    // Choose a reference tile
    int refinp = 84; // Tile012 (?)
    Antenna ref_ant = obs_metadata->antennas[refinp];
    double N_ref = ref_ant.north_m;
    double E_ref = ref_ant.east_m;
    double H_ref = ref_ant.height_m;

    long int freq_ch;

    double cable;

    double integer_phase;
    double geometry, delay_time, delay_samples, cycles_per_sample;

    double uv_angle;
    cuDoubleComplex uv_phase; // For the UV phase correction

    double Fnorm;

    for (int p = 0; p < npointing; p++)
    {
        // Everything from this point on is frequency-dependent
        for (ch = 0; ch < nchan; ch++) {

            // Calculating direction-dependent matrices
            freq_ch = frequency + ch*chan_width;    // The frequency of this fine channel

            // Calculate the UV phase correction for this channel
            uv_angle = cal->phase_slope*freq_ch + cal->phase_offset;
            uv_phase = make_cuDoubleComplex( cos(uv_angle), sin(uv_angle) );

            for (rf_input = 0; rf_input < (int)(ninput); rf_input++) {

                // Get the antenna and polarisation number from the rf_input
                ant = obs_metadata->rf_inputs[rf_input].ant;
                pol = *(obs_metadata->rf_inputs[rf_input].pol) - 'X'; // 'X' --> 0; 'Y' --> 1

                mult2x2d(D[ant][ch], B[p][ant], Ji); // the gain in the desired look direction

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
    } // end loop through pointings (p)

}

void calc_beam_geom(
        double           *ras_hours,
        double           *decs_degs,
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

    double pr = 0, pd = 0, px = 0, rv = 0, eq = 2000, ra_ap = 0, dec_ap = 0;

    /* get mjd */
    intmjd = floor(mjd);
    fracmjd = mjd - intmjd;

    /* get requested Az/El from command line */
    mjd2lst( mjd, &lmst );

    for (int p = 0; p < npointing; p++)
    {
        /* for the look direction <not the tile> */

        mean_ra = ras_hours[p] * PAL__DH2R;
        mean_dec = decs_degs[p] * PAL__DD2R;

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
