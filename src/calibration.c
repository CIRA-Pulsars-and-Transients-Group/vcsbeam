/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <mwalib.h>
#include <cuComplex.h>

#include "vcsbeam.h"

cuDoubleComplex *get_rts_solution( MetafitsMetadata *cal_metadata,
        MetafitsMetadata *obs_metadata, const char *caldir, uintptr_t coarse_chan_idx )
/* Read in the RTS solution from the DI_Jones... and Bandpass... files in
 * the CALDIR directory. The output is a set of Jones matrices (D) for each
 * antenna and (non-flagged) fine channel.
 *
 * This function allocates memory, which should be freed by the caller (free())
 */
{
    // Find the "GPUBox" number for this coarse channel
    uintptr_t gpubox_number = cal_metadata->metafits_coarse_chans[coarse_chan_idx].gpubox_number;

    // With the gpubox number in hand, construct the filenames for the
    // DI_Jones and Bandpass files
    char dijones_path[CAL_BUFSIZE];
    char bandpass_path[CAL_BUFSIZE];

    sprintf( dijones_path,  "%s/DI_JonesMatrices_node%03lu.dat", caldir, gpubox_number );
    sprintf( bandpass_path, "%s/BandpassCalibration_node%03lu.dat", caldir, gpubox_number );

    // Allocate memory for the Jones arrays
    uintptr_t nant    = cal_metadata->num_ants;
    uintptr_t nvispol = cal_metadata->num_visibility_pols; // = 4 (XX, XY, YX, YY)
    uintptr_t nantpol = cal_metadata->num_ant_pols; // = 2 (X, Y)
    uintptr_t nchan   = cal_metadata->num_corr_fine_chans_per_coarse;

    // Temporary arrays for DI Jones matrices ('Dd') and Bandpass matrices ('Db')
    cuDoubleComplex  **Dd = (cuDoubleComplex ** )calloc( nant, sizeof(cuDoubleComplex * ) );
    cuDoubleComplex ***Db = (cuDoubleComplex ***)calloc( nant, sizeof(cuDoubleComplex **) );

    // The array for the final output product (D = Dd x Db)
    // This will have the same dimensions as the final Jones matrix, so use
    // J_IDX for indexing
    size_t Dsize = nant*nchan*nvispol;
    cuDoubleComplex *D  = (cuDoubleComplex *)calloc( Dsize, sizeof(cuDoubleComplex) );

    uintptr_t ant, ch;
    for (ant = 0; ant < nant; ant++)
    {
        Dd[ant] = (cuDoubleComplex * )calloc( nvispol, sizeof(cuDoubleComplex  ) );
        Db[ant] = (cuDoubleComplex **)calloc( nchan,   sizeof(cuDoubleComplex *) );

        for (ch = 0; ch < nchan; ch++)
        {
            Db[ant][ch] = (cuDoubleComplex *)calloc( nvispol, sizeof(cuDoubleComplex) );
        }
    }

    // Read in the DI Jones file
    read_dijones_file((double **)Dd, NULL, cal_metadata->num_ants, dijones_path);

    // Read in the Bandpass file
    read_bandpass_file( NULL, Db, cal_metadata, bandpass_path );

    // Form the "fine channel" DI gain (the "D" in Eqs. (28-30), Ord et al. (2019))
    // Do this for each of the _voltage_ observation's fine channels (use
    // nearest-neighbour interpolation). This is "applying the bandpass" corrections.
    uintptr_t vcs_nchan = obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t interp_factor = vcs_nchan / nchan;
    uintptr_t cal_ch, cal_ant;
    uintptr_t d_idx;
    for (ant = 0; ant < nant; ant++)
    {
        // Match up antenna indices between the cal and obs metadata
        cal_ant = get_idx_for_vcs_antenna_in_cal( cal_metadata, obs_metadata, ant );

        for (ch = 0; ch < vcs_nchan; ch++)
        {
            // A pointer to the (first element of) the Jones matrix for this
            // antenna and channel (d_idx) and the Jones matrix for the first
            // antenna (our reference antenna)
            d_idx    = J_IDX(ant,ch,0,0,nchan,nantpol);

            // D = Dd x Db
            // SM: Daniel Mitchell confirmed in an email dated 23 Mar 2017 that the
            // Bandpass matrices (Db) should be multiplied on the _right_ of the
            // DI Jones matrices (Dd).
            cal_ch = ch / interp_factor;
            mult2x2d( Dd[cal_ant], Db[cal_ant][cal_ch], &(D[d_idx]) );
            cp2x2( Dd[cal_ant], &(D[d_idx]) );
        }
    }

    // Free memory for temporary arrays
    for (ant = 0; ant < nant; ant++)
    {
        for (ch = 0; ch < nchan; ch++)
            free( Db[ant][ch] );
        free( Db[ant] );
        free( Dd[ant] );
    }

    free( Db );
    free( Dd );

    return D;
}

void read_dijones_file( double **Dd, double *amp, uintptr_t nant, char *fname )
/* Read in an RTS file and return the direction independent Jones matrix
 * for each antenna. This implements Eq. (29) in Ord et al. (2019).
 *
 * This function assumes that the RTS "DIJones" files are in a very specific
 * format. However, this code is maintained independently from the RTS,
 * so if the RTS changes, this code may break.
 */
{
    // Open the file for reading
    FILE *fp = NULL;
    if ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open gain Jones matrix file: %s\n",
                fname);
        exit(EXIT_FAILURE);
    }

    // Set up some variables for reading in
    double val;
    double invA[NDBL_PER_JONES]; // The alignment matrix
    double J[NDBL_PER_JONES]; // "Raw" Jones matrix

    // Read in the amplitude (SM: I have no idea what this value represents)
    fscanf( fp, "%lf", &val );
    if (amp != NULL)
        *amp = val;

    // Read in the alignment ("A") matrix and invert it ("invA")
    // (Just use the "invA" array for both reading and inverting)

    // Reading:
    int i;
    for (i = 0; i < NDBL_PER_JONES; i++)
        fscanf( fp, "%lf,", &invA[i] );

    // Inverting (inv2x2() expects cuDoubleComplex arrays):
    inv2x2( (cuDoubleComplex *)invA, (cuDoubleComplex *)invA );

    // Finally, read in the Jones ("J") matrices and multiply them each by the
    // inverted alignment matrix ("invA")
    // Eq. (29) in Ord et al. (2019)
    uintptr_t ant;
    for (ant = 0; ant < nant; ant++)
    {
        for (i = 0; i < NDBL_PER_JONES; i++)
            fscanf( fp, "%lf,", &J[i] );

        mult2x2d( (cuDoubleComplex *)J, (cuDoubleComplex *)invA,
                (cuDoubleComplex *)(Dd[ant]) );
    }

    fclose(fp);
}



void read_bandpass_file(
        cuDoubleComplex ***Jm, // Output: measured Jones matrices (Jm[ant][ch][pol,pol])
        cuDoubleComplex ***Jf, // Output: fitted Jones matrices   (Jf[ant][ch][pol,pol])
        MetafitsMetadata  *cal_metadata,
        char *filename         // Input:  name of bandpass file
        )
/* This function populates the Jm and Jf arrays with values read in from the
 * given bandpass file. The bandpass files contain only values for antennas
 * and fine channels that have not been flagged. Nothing is done for those
 * antennas/channels that have been flagged, so the onus is on the caller to
 * initialise the Jm and Jf arrays to values to their preferred values.
 */
{

    // Some shortcut variables
    int chan_width = cal_metadata->corr_fine_chan_width_hz;        // TODO: Check whether this should sometimes be volt_fine_chan_width_hz
    int nchan      = cal_metadata->num_corr_fine_chans_per_coarse; //       and num_volt_fine_chans_per_coarse
    int nant       = cal_metadata->num_ants;

    // Open the file for reading
    FILE *f = NULL;
    if ((f = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open Bandpass file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Read in top row = frequency offsets within coarse channel
    // Turn these frequencies into a list of indexes into the channel dimension
    char freqline[CAL_BUFSIZE];
    if (fgets( freqline, CAL_BUFSIZE, f ) == NULL) {
        fprintf(stderr, "Error: could not read first line of %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Parse top row
    // Find out which channels are actually in the Bandpass file
    // (i.e. which channels have not been flagged)
    int chan_idxs[nchan];
    char *freqline_ptr = freqline;
    int pos;
    double chan_MHz;
    int chan_count = 0;
    int chan_idx;
    while (sscanf( freqline_ptr, "%lf,%n", &chan_MHz, &pos ) == 1)
    {
        chan_count++;

        // Make sure we haven't exceeded the total
        if (chan_count > nchan) {
            fprintf(stderr, "Error: More than nchan = %d columns in Bandpass file %s\n", nchan, filename);
            exit(EXIT_FAILURE);
        }

        // Convert the channel frequencies (in MHz) to channel indices
        chan_idx = (int)roundf( chan_MHz*1e6 / (double)chan_width );
        chan_idxs[chan_count - 1] = chan_idx;

        freqline_ptr += pos;
    }

    // Read in all the values
    int ant, curr_ant;         // The antenna number, read in from the file
    int ant_row;               // Number between 0 and 7. Each antenna has 8 rows.
    int ch;                    // A counter for channel numbers
    int ci;                    // A counter for channel number indices
    double amp,ph;             // For holding the read-in value pairs (real, imaginary)
    int pol;                   // Number between 0 and 3. Corresponds to position in Jm/Jf matrices: [0,1]
                               //                                                                    [2,3]
    cuDoubleComplex ***J;      // Either points to Jm or Jf, according to which row we're on

    for (ant = 0; ant < nant; ant++)
    {
        // Read in first number = antenna number + 1
        if (fscanf( f, "%d,", &curr_ant ) == EOF)
            break;

        // (The antenna numbers in the Bandpass files start at 1, not 0)
        curr_ant--;

        // If the current row is not for this ant, then skip to the current row
        if (ant < curr_ant)  ant = curr_ant;

        // Now we've caught up to the next non-absent antenna in the bandpass file, so
        // go through the next 8 rows, which should all pertain to this antenna
        for (ant_row = 0; ant_row < BANDPASS_ROWS_PER_ANT; ant_row++)
        {
            // If we're on the first of 8 lines, then we've already consumed
            // the antenna number at the beginning of that row. Otherwise, we
            // still need to do it. We will re-purpose "curr_ant" for this.
            if (ant_row > 0)
            {
                fscanf( f, "%d,", &curr_ant );
                curr_ant--;

                // Sanity check: the antenna number hasn't changed within the 8 rows
                if (ant != curr_ant)
                {
                    fprintf( stderr, "error: read_bandpass_file: fewer than eight rows "
                            "detected for antenna %u (labelled '%u') in %s\n",
                            ant, ant + 1, filename );
                    exit(EXIT_FAILURE);
                }
            }

            // Decide if the row corresponds to the Jm values (even rows)
            // or Jf values (odd rows)
            if (ant_row % 2 == 0)
                J = Jm;
            else
                J = Jf;

            // If the caller doesn't care about this row, skip the rest of
            // this line (freqline isn't needed any more, so use it as a dummy
            // buffer) and start afresh on the next line
            if (J == NULL)
            {
                fgets( freqline, CAL_BUFSIZE, f );
                continue;
            }

            // Get the polarisation index
            pol = ant_row/2;

            // Loop over the row and read in values
            for (ci = 0; ci < chan_count; ci++)
            {
                ch = chan_idxs[ci];                 // Get the channel number
                fscanf( f, "%lf,%lf,", &amp, &ph ); // Read in the re,im pairs in each row

                // Convert to complex number and store in output array
                J[ant][ch][pol] = make_cuDoubleComplex( amp*cos(ph), amp*sin(ph) );
            }
        } // End for loop over 8 rows (for one antenna)
    }
}


int read_offringa_gains_file( cuDoubleComplex **antenna_gain, int nant,
                              int coarse_chan, char *gains_file )
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
    int pol_idx, ant_idx;   // Used for "re-ordering" the antennas and pols (<-- TODO: FIX ME)
    int count = 0;          // Keep track of how many solutions have actually been read in
    double re, im;          // Temporary placeholders for the real and imaginary doubles read in

    // Loop through antennas and read in calibration solution
    int first = 1;
    for (ant = 0; ant < nant; ant++) {

        ant_idx = ant;

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


void remove_reference_phase( cuDoubleComplex *J, cuDoubleComplex *Jref )
{
    cuDoubleComplex XX0norm, XY0norm, YX0norm, YY0norm;
    double XXscale = 1.0/cuCabs( Jref[0] ); // = 1/|XX|
    double XYscale = 1.0/cuCabs( Jref[1] ); // = 1/|XY|
    double YXscale = 1.0/cuCabs( Jref[2] ); // = 1/|YX|
    double YYscale = 1.0/cuCabs( Jref[3] ); // = 1/|YY|

    XX0norm = make_cuDoubleComplex( XXscale*cuCreal(Jref[0]), XXscale*cuCimag(Jref[0]) ); // = XX/|XX|
    XY0norm = make_cuDoubleComplex( XYscale*cuCreal(Jref[1]), XYscale*cuCimag(Jref[1]) ); // = XY/|XY|
    YX0norm = make_cuDoubleComplex( YXscale*cuCreal(Jref[2]), YXscale*cuCimag(Jref[2]) ); // = YX/|YX|
    YY0norm = make_cuDoubleComplex( YYscale*cuCreal(Jref[3]), YYscale*cuCimag(Jref[3]) ); // = YY/|YY|

    J[0] = cuCdiv( J[0], XX0norm ); // Essentially phase rotations
    J[1] = cuCdiv( J[1], XY0norm );
    J[2] = cuCdiv( J[2], YX0norm );
    J[3] = cuCdiv( J[3], YY0norm );
}

void zero_XY_and_YX( cuDoubleComplex *J )
/* For J = [ XX, XY ], set XY and YX to 0 for all antennas
 *         [ YX, YY ]
 */
{
    J[1] = make_cuDoubleComplex( 0.0, 0.0 );
    J[2] = make_cuDoubleComplex( 0.0, 0.0 );
}

void pq_phase_correction( uint32_t gpstime, cuDoubleComplex *D, MetafitsMetadata *obs_metadata, logger *log )
/* Retrieve the XY phase correction from pq_phase_correction.txt for the given
 * gps time
 */
{
    char log_message[CAL_BUFSIZE + 100];
    char buffer[CAL_BUFSIZE];
    sprintf( buffer, "%s/pq_phase_correction.txt", RUNTIME_DIR );
    FILE *f = fopen( buffer, "r" );
    if (f == NULL)
    {
        sprintf( log_message, "warning: Cannot find file '%s' -- will not apply any PQ correction",
                buffer );
        logger_timed_message( log, log_message );

        return;
    }

    // Values to be read in
    uint32_t from, to;
    double slope = 0.0, offset = 0.0;
    char ref_tile_name[32];

    // Read in the file line by line
    while (fgets( buffer, CAL_BUFSIZE, f ) != NULL)
    {
        // Only parse if the first character is numeric
        if (!(buffer[0] >= '0' && buffer[0] <= '9'))
            continue;

        // Parse the line
        if (sscanf( buffer, "%u %u %lf %lf %s", &from, &to, &slope, &offset, ref_tile_name ) != 5)
        {
            logger_timed_message( log, "ERROR: pq_phase_correction: cannot parse PQ phase correction file" );
            exit(EXIT_FAILURE);
        }

        // If the given gpstime is in the right range, keep the values and
        // stop looking
        if (from <= gpstime && gpstime <= to)
            break;
    }

    // Close the file
    fclose( f );

    if (slope == 0.0 && offset == 0.0)
    {
        // Nothing else to do, so exit
        return;
    }

    // Get the reference antenna index
    int ref_ant = find_antenna_by_name( obs_metadata, ref_tile_name );
    if (ref_ant == NO_ANTENNA_FOUND)
    {
        sprintf( log_message, "PQ-PHASE: Antenna \"%s\" not found in this observation", ref_tile_name );
        logger_timed_message( log, log_message );

        ref_ant = 0;
        sprintf( log_message, "PQ-PHASE: Using (default) \"%s\" (antenna #0 in VCS) instead",
                obs_metadata->antennas[ref_ant].tile_name );
        logger_timed_message( log, log_message );
    }
    uintptr_t d_idx, dref_idx;

    uintptr_t nant    = obs_metadata->num_ants;
    uintptr_t nantpol = obs_metadata->num_ant_pols; // = 2 (X, Y)
    uintptr_t nchan   = obs_metadata->num_corr_fine_chans_per_coarse;

    // A temporary copy of the reference antenna matrix, so that we don't clobber it midway
    // through this operation!
    cuDoubleComplex Dref[nantpol*nantpol];

    // Go through all the antennas and divide the phases by the reference antenna
    uintptr_t ant, ch;
    for (ch = 0; ch < nchan; ch++)
    {
        // Make a copy of the reference Jones matrix for this channel
        // reference antenna and channel
        dref_idx = J_IDX(ref_ant,ch,0,0,nchan,nantpol);
        cp2x2( &(D[dref_idx]), Dref );

        for (ant = 0; ant < nant; ant++)
        {
            // A pointer to the (first element of) the Jones matrix for this
            // antenna and channel
            d_idx = J_IDX(ant,ch,0,0,nchan,nantpol);

            // By default, divide through a reference antenna...
            remove_reference_phase( &(D[d_idx]), Dref );

            // ...and zero the off-diagonal terms
            zero_XY_and_YX( &(D[d_idx]) );
        }
    }
}
