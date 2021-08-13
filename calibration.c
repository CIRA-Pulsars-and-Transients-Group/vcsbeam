/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include "calibration.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mwalib.h>
#include <cuComplex.h>

int read_rts_file(cuDoubleComplex **G, cuDoubleComplex *Jref,
                  double *amp, char *fname)
{
    FILE *fp = NULL;
    if ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open gain Jones matrix file: %s\n",
                fname);
        exit(EXIT_FAILURE);
    }

    char line[BUFSIZE];
    int index = 0;
    double re0, im0, re1, im1, re2, im2, re3, im3;

    while ((fgets(line, BUFSIZE - 1, fp)) != NULL) {

        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0')
            continue; // skip blank/comment lines
        if (line[0] == '/' && line[1] == '/')
            continue; // also a comment (to match other input files using this style)

        if (index == 0) {

            // read the amplitude and the Alignment Line
            sscanf(line, "%lf", amp);
            fgets(line, BUFSIZE - 1, fp);
            sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);

            Jref[0] = make_cuDoubleComplex( re0, im0 );
            Jref[1] = make_cuDoubleComplex( re1, im1 );
            Jref[2] = make_cuDoubleComplex( re2, im2 );
            Jref[3] = make_cuDoubleComplex( re3, im3 );

        }
        if (index > 0) {
            sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);
            G[index - 1][0] = make_cuDoubleComplex( re0, im0 );
            G[index - 1][1] = make_cuDoubleComplex( re1, im1 );
            G[index - 1][2] = make_cuDoubleComplex( re2, im2 );
            G[index - 1][3] = make_cuDoubleComplex( re3, im3 );
        }

        index++;

    }

    fclose(fp);

    return 0;

}



int read_bandpass_file(
        cuDoubleComplex ***Jm, // Output: measured Jones matrices (Jm[ant][ch][pol,pol])
        cuDoubleComplex ***Jf, // Output: fitted Jones matrices   (Jf[ant][ch][pol,pol])
        int chan_width,       // Input:  channel width of one column in file (in Hz)
        int nchan,            // Input:  (max) number of channels in one file (=128/(chan_width/10000))
        int nant,             // Input:  (max) number of antennas in one file (=128)
        char *filename        // Input:  name of bandpass file
        )
{

    // Open the file for reading
    FILE *f = NULL;
    if ((f = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open Bandpass file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Read in top row = frequency offsets
    int max_len = 4096; // Overkill
    char freqline[max_len];
    if (fgets( freqline, max_len, f ) == NULL) {
        fprintf(stderr, "Error: could not read first line of %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Parse top row
    // Find out which channels are actually in the Bandpass file
    // (i.e. which channels have not been flagged)
    char *freqline_ptr = freqline;
    int pos;
    double freq_offset;
    int chan_count = 0;
    int chan_idxs[nchan];
    int chan_idx;
    while (sscanf(freqline_ptr, "%lf,%n", &freq_offset, &pos) == 1) {

        chan_count++;

        // Make sure we haven't exceeded the total
        if (chan_count > nchan) {
            fprintf(stderr, "Error: More than nchan = %d columns in Bandpass file %s\n", nchan, filename);
            exit(EXIT_FAILURE);
        }
        chan_idx = (int)roundf( freq_offset*1e6 / (double)chan_width );
        chan_idxs[chan_count-1] = chan_idx;

        freqline_ptr += pos;
    }
    // Read in all the values
    int ant, curr_ant = 0;     // The antenna number, read in from the file
    int ant_row       = 0;     // Number between 0 and 7. Each antenna has 8 rows.
    int ch;                    // A counter for channel numbers
    int ci;                    // A counter for channel number indices
    double amp,ph;             // For holding the read-in value pairs (real, imaginary)
    int pol;                   // Number between 0 and 3. Corresponds to position in Jm/Jf matrices: [0,1]
                               //                                                                    [2,3]
    cuDoubleComplex ***J;       // Either points to Jm or Jf, according to which row we're on

    while (1) {   // Will terminate when EOF is reached

        if (fscanf(f, "%d,", &ant) == EOF)     // Read in first number = antenna number
            break;

        if (ant > nant) {                      // Check that the antenna number is not bigger than expected
            fprintf(stderr, "Error: More than nant = %d antennas in Bandpass file %s\n", nant, filename);
            exit(EXIT_FAILURE);
        }

        ant--;                                 // Convert to 0-offset

        if (ant == curr_ant) {                 // Ensure that there is not an unusual (!=8) number of rows for this antenna
            ant_row++;
            if (ant_row > 8) {
                fprintf(stderr, "Error: More than 8 rows for antenna %d in Bandpass file %s\n",
                        ant, filename);
                exit(EXIT_FAILURE);
            }
        }
        else {
            if (ant_row < 7) {
                fprintf(stderr, "Error: Fewer than 8 rows for antenna %d in Bandpass file %s\n",
                        ant, filename);
                exit(EXIT_FAILURE);
            }
            curr_ant = ant;
            ant_row  = 1;
        }

        if ((ant_row-1) % 2 == 0)  J = Jm;      // Decide if the row corresponds to the Jm values (even rows)
        else                       J = Jf;      // or Jf values (odd rows)

        if (J == NULL) {                        // If the caller doesn't care about this row
            fgets( freqline, max_len, f );      // Skip the rest of this line (freqline isn't needed any more)
            continue;                           // And start afresh on the next line
        }

        pol = (ant_row-1) / 2;                  // Get the polarisation index

        for (ci = 0; ci < chan_count; ci++) {   // Loop over the row

            ch = chan_idxs[ci];                 // Get the channel number
            fscanf(f, "%lf,%lf,", &amp, &ph);   // Read in the re,im pairs in each row

            J[ant][ch][pol] = make_cuDoubleComplex( amp*cos(ph), amp*sin(ph) );
                                                // Convert to complex number and store in output array
        }

        // (Assumes that the number of values in each row are correct)
    }

    return 1;
}


int read_offringa_gains_file( cuDoubleComplex **antenna_gain, int nant,
                              int coarse_chan, char *gains_file, int *order )
{
    // Assumes that memory for antenna has already been allocated

    // Open the calibration file for reading
    FILE *fp = NULL;
    fp = fopen(gains_file,"r");
    if (fp == NULL) {
        fprintf(stderr,"Failed to open %s: quitting\n",gains_file);
        exit(EXIT_FAILURE);
    }

    // Read in the necessary information from the header

    uint32_t intervalCount, antennaCount, channelCount, polarizationCount;

    fseek(fp, 16, SEEK_SET);
    fread(&intervalCount,     sizeof(uint32_t), 1, fp);
    fread(&antennaCount,      sizeof(uint32_t), 1, fp);
    fread(&channelCount,      sizeof(uint32_t), 1, fp);
    fread(&polarizationCount, sizeof(uint32_t), 1, fp);

    // Error-checking the info extracted from the header
    if (intervalCount > 1) {
        fprintf(stderr, "Warning: Only the first interval in the calibration ");
        fprintf(stderr, "solution (%s) will be used\n", gains_file);
    }
    if ((int)antennaCount != nant) {
        fprintf(stderr, "Error: Calibration solution (%s) ", gains_file);
        fprintf(stderr, "contains a different number of antennas (%d) ", antennaCount);
        fprintf(stderr, "than specified (%d)\n", nant);
        exit(1);
    }
    if (channelCount != 24) {
        fprintf(stderr, "Warning: Calibration solution (%s) ", gains_file);
        fprintf(stderr, "contains a different number (%d) ", channelCount);
        fprintf(stderr, "than the expected (%d) channels. ", 24);
    }
    if ((int)channelCount <= coarse_chan) {
        fprintf(stderr, "Error: Requested channel number (%d) ", coarse_chan);
        fprintf(stderr, "is more than the number of channels (0-%d) ", channelCount-1);
        fprintf(stderr, "available in the calibration solution (%s)\n", gains_file);
        exit(1);
    }
    int npols = polarizationCount; // This will always = 4

    // Prepare to jump to the first solution to be read in
    int bytes_left_in_header = 16;
    int bytes_to_first_jones = bytes_left_in_header + (npols * coarse_chan * sizeof(cuDoubleComplex));
         //     (See Offringa's specs for details)
         //     Assumes coarse_chan is zero-offset
         //     sizeof(complex double) *must* be 64-bit x 2 = 16-byte
    int bytes_to_next_jones = npols * (channelCount-1) * sizeof(cuDoubleComplex);

    int ant, pol;           // Iterate through antennas and polarisations
    int pol_idx, ant_idx;   // Used for "re-ordering" the antennas and pols
    int count = 0;          // Keep track of how many solutions have actually been read in
    double re, im;          // Temporary placeholders for the real and imaginary doubles read in

    // Loop through antennas and read in calibration solution
    int first = 1;
    for (ant = 0; ant < nant; ant++) {

        // Get correct antenna index
        // To wit: The nth antenna in the Offringa binary file will get put into
        // position number order[n]. Default is no re-ordering.
        if (order) {
            ant_idx = order[ant];
        }
        else {
            ant_idx = ant;
        }

        // Jump to next Jones matrix position for this channel
        if (first) {
            fseek(fp, bytes_to_first_jones, SEEK_CUR);
            first = 0;
        }
        else {
            fseek(fp, bytes_to_next_jones, SEEK_CUR);
        }

        // Read in the data
        for (pol = 0; pol < npols; pol++) {

            pol_idx = 3 - pol; // Read them in "backwards", because RTS's "x" = Offringa's "y"

            fread(&re, sizeof(double), 1, fp);
            fread(&im, sizeof(double), 1, fp);

            // Check for NaNs
            if (isnan(re) | isnan(im)) {

                // If NaN, set to identity matrix
                if (pol_idx == 0 || pol_idx == 3)
                    antenna_gain[ant_idx][pol_idx] = make_cuDoubleComplex( 1.0, 0.0 );
                else
                    antenna_gain[ant_idx][pol_idx] = make_cuDoubleComplex( 0.0, 0.0 );

            }
            else {
                antenna_gain[ant_idx][pol_idx] = make_cuDoubleComplex( re, im );
            }

            count++;

        }
    }

    // Close the file, print a summary, and return
    fclose(fp);
    fprintf(stdout, "Read %d inputs from %s\n", count, gains_file);

    return count/npols; // should equal the number of antennas
}


uint32_t get_idx_for_vcs_antenna_in_cal( MetafitsMetadata *cal_metadata, MetafitsMetadata *obs_metadata, uint32_t vcs_ant )
/* Get the (antenna) idx into the calibration solutions that corresponds to the supplied vcs ANTenna number.
 * Because this is attempting to match up antennas across two different observation, there is no
 * guarantee that the two observations will even use the same sets of antennas (i.e. tiles). Thus,
 * in order for this function to be robust, antennas must be matched by "tile_id", which is uniquely
 * assigned to each physical tile.
 */
{
    // Find the tile_id for this antenna number in the vcs observation ("obs")
    uint32_t tile_id, i;
    bool found = false;
    for (i = 0; i < obs_metadata->num_ants; i++)
    {
        if (obs_metadata->antennas[i].ant == vcs_ant)
        {
            tile_id = obs_metadata->antennas[i].tile_id;
            found = true;
            break;
        }
    }

    // If they've provided a non-existent vcs_ant, then complain and exit
    if (!found)
    {
        fprintf( stderr, "error: get_idx_for_vcs_antenna_in_cal: could not "
                "find antenna %u in vcs obs %u\n",
                vcs_ant, obs_metadata->obs_id );
        exit(EXIT_FAILURE);
    }

    // Now find the matching idx for this tile_id in the calibration observation ("cal")
    uint32_t idx;
    found = false;
    for (idx = 0; idx < cal_metadata->num_ants; idx++)
    {
        if (cal_metadata->antennas[idx].tile_id == tile_id)
        {
            found = true;
            break;
        }
    }

    // If it was not found at all, then exit with an error
    if (!found)
    {
        fprintf( stderr, "error: get_idx_for_vcs_antenna_in_cal: could not "
                "find tile_id %u in cal obs %u, required by vcs obs %u\n",
                tile_id, cal_metadata->obs_id, obs_metadata->obs_id );
        exit(EXIT_FAILURE);
    }

    // Otherwise, all good!
    return cal_metadata->rf_inputs[cal_metadata->antennas[idx].rfinput_x].input/2;
}
