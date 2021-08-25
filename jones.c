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
#include "jones.h"
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

void get_jones(
        // an array of pointings [pointing][ra/dec][characters]
        int                    npointing, // number of pointings
        MetafitsMetadata      *obs_metadata,
        int                    coarse_chan_idx,
        struct                 calibration *cal,
        cuDoubleComplex     ***D,
        cuDoubleComplex       *B,
        cuDoubleComplex    ****invJi )                 // output: invJi[ant][ch][pol][pol]
{

    // Give "shorthand" variables for often-used values in metafits
    long int frequency = obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_start_hz;
    int nant           = obs_metadata->num_ants;
    int nchan          = obs_metadata->num_volt_fine_chans_per_coarse;
    int npol           = obs_metadata->num_ant_pols;   // (X,Y)
    int chan_width     = obs_metadata->corr_fine_chan_width_hz;

    int ant;     // Antenna number
    int ch;      // Channel number
    int p1, p2;  // Counters for polarisation

    /* easy -- now the positions from the database */

    cuDoubleComplex Ji[npol*npol];              // Gain in Desired Direction

    long int freq_ch;

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

            for (ant = 0; ant < nant; ant++)
            {
                mult2x2d(D[ant][ch], &(B[PB_IDX(p, ant, 0, nant, npol*npol)]), Ji); // the gain in the desired look direction

                // Apply the UV phase correction (to the bottom row of the Jones matrix)
                Ji[2] = cuCmul( Ji[2], uv_phase );
                Ji[3] = cuCmul( Ji[3], uv_phase );

                // Now, calculate the inverse Jones matrix

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

            } // end loop through antenna/pol (rf_input)
        } // end loop through fine channels (ch)
    } // end loop through pointings (p)

}

