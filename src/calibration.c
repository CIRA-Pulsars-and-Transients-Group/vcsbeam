/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <mwalib.h>
#include <cuComplex.h>

#include "vcsbeam.h"

cuDoubleComplex *get_rts_solution( MetafitsMetadata *cal_metadata,
        MetafitsMetadata *obs_metadata, const char *caldir, uintptr_t coarse_chan_idx, logger *log )
/* Read in the RTS solution from the DI_Jones... and Bandpass... files in
 * the CALDIR directory. The output is a set of Jones matrices (D) for each
 * antenna and (non-flagged) fine channel.
 *
 * This function allocates memory, which should be freed by the caller (free())
 */
{
    // Find the "GPUBox" number for this coarse channel
    //uintptr_t gpubox_number = cal_metadata->metafits_coarse_chans[coarse_chan_idx].gpubox_number; // <-- This is what it should be, when the mwalib bug is fixed
    uintptr_t gpubox_number = cal_metadata->metafits_coarse_chans[coarse_chan_idx].corr_chan_number + 1; // <-- This is the temporary hack

    // With the gpubox number in hand, construct the filenames for the
    // DI_Jones and Bandpass files
    char dijones_path[CAL_BUFSIZE];
    char bandpass_path[CAL_BUFSIZE];

    sprintf( dijones_path,  "%s/DI_JonesMatrices_node%03lu.dat", caldir, gpubox_number );
    sprintf( bandpass_path, "%s/BandpassCalibration_node%03lu.dat", caldir, gpubox_number );

    // Allocate memory for the Jones arrays
    uintptr_t nant    = cal_metadata->num_ants;
    uintptr_t ninputs = cal_metadata->num_rf_inputs;
    uintptr_t nvispol = cal_metadata->num_visibility_pols; // = 4 (PP, PQ, QP, QQ)
    uintptr_t nantpol = cal_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t nchan   = cal_metadata->num_corr_fine_chans_per_coarse;
    uintptr_t vcs_nchan = obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t interp_factor = vcs_nchan / nchan;
    uintptr_t ant, ch; // For loop variables

    // Write out the channel --> gpubox number mapping
    char log_message[256];
    if (log)
    {
        sprintf( log_message, "Receiver channel #%lu --> GPUBox #%lu",
                cal_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number,
                gpubox_number );
        logger_timed_message( log, log_message );
    }

    // Temporary arrays for DI Jones matrices ('Dd') and Bandpass matrices ('Db')
    cuDoubleComplex  **Dd = (cuDoubleComplex ** )calloc( nant, sizeof(cuDoubleComplex * ) );
    cuDoubleComplex ***Db = (cuDoubleComplex ***)calloc( nant, sizeof(cuDoubleComplex **) );

    // The array for the final output product (D = Dd x Db)
    // This will have the same dimensions as the final Jones matrix, so use
    // J_IDX for indexing
    size_t Dsize = nant*vcs_nchan*nvispol;
    cuDoubleComplex *D  = (cuDoubleComplex *)calloc( Dsize, sizeof(cuDoubleComplex) );

    for (ant = 0; ant < nant; ant++)
    {
        Dd[ant] = (cuDoubleComplex * )calloc( nvispol, sizeof(cuDoubleComplex  ) );
        Db[ant] = (cuDoubleComplex **)calloc( nchan,   sizeof(cuDoubleComplex *) );

        for (ch = 0; ch < nchan; ch++)
            Db[ant][ch] = (cuDoubleComplex *)calloc( nvispol, sizeof(cuDoubleComplex) );
    }

    // Read in the DI Jones file
    read_dijones_file((double **)Dd, NULL, cal_metadata->num_ants, dijones_path);

    // Read in the Bandpass file
    read_bandpass_file( NULL, Db, cal_metadata, bandpass_path );

    // Make the master mpi thread print out the antenna names of both
    // obs and cal metafits. "Header" printed here, actual numbers
    // printed inside antenna for loop below
    if (log)
    {
        if (log->world_rank == 0)
        {
            logger_message( log, "" );
            logger_timed_message( log, "Calibration metafits info:" );
            logger_timed_message( log, "--------------------------" );
            logger_timed_message( log, "|Input| Ant |Tile ID| TileName |Pol|VCS order|" );

            uintptr_t i;
            for (i = 0; i < ninputs; i++)
            {
                sprintf( log_message, "| %3u | %3u | %5u | %8s | %c |  %3u   |",
                        cal_metadata->rf_inputs[i].input,
                        cal_metadata->rf_inputs[i].ant,
                        cal_metadata->rf_inputs[i].tile_id,
                        cal_metadata->rf_inputs[i].tile_name,
                        *(cal_metadata->rf_inputs[i].pol),
                        cal_metadata->rf_inputs[i].vcs_order
                       );
                logger_timed_message( log, log_message );
            }

            logger_message( log, "" );
        }
    }

    // Form the "fine channel" DI gain (the "D" in Eqs. (28-30), Ord et al. (2019))
    // Do this for each of the _voltage_ observation's fine channels (use
    // nearest-neighbour interpolation). This is "applying the bandpass" corrections.
    uintptr_t cal_ch, cal_ant, obs_ant, dd_idx;
    uintptr_t d_idx;
    uintptr_t i; // (i)dx into rf_inputs
    Rfinput *obs_rfinput, *cal_rfinput;
    for (i = 0; i < ninputs; i++)
    {
        // Order the Jones matrices by the TARGET OBSERVATION metadata
        // (not the calibration metadata)
        obs_rfinput = &(obs_metadata->rf_inputs[i]);

        // Only have to loop once per tile, so skip the 'Y's
        if (*(obs_rfinput->pol) == 'Y')
            continue;

        // Match up antenna indices between the cal and obs metadata
        cal_rfinput = find_matching_rf_input( cal_metadata, obs_rfinput );

        obs_ant = obs_rfinput->ant;
        dd_idx  = cal_rfinput->input/2;

        for (ch = 0; ch < vcs_nchan; ch++)
        {
            // A pointer to the (first element of) the Jones matrix for this
            // antenna and channel (d_idx) and the Jones matrix for the first
            // antenna (our reference antenna)
            d_idx = J_IDX(obs_ant,ch,0,0,vcs_nchan,nantpol);

            // D = Dd x Db
            // SM: Daniel Mitchell confirmed in an email dated 23 Mar 2017 that the
            // Bandpass matrices (Db) should be multiplied on the _right_ of the
            // DI Jones matrices (Dd).
            cal_ch = ch / interp_factor;
            mult2x2d( Dd[dd_idx], Db[dd_idx][cal_ch], &(D[d_idx]) );
            //cp2x2( Dd[dd_idx], &(D[d_idx]) );
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
//fprintf( stderr, "chan_MHz   = %lf\n", chan_MHz );
//fprintf( stderr, "chan_width = %lf\n", (double)chan_width/1e6 );
//fprintf( stderr, "chan_idx   = %d\n\n", chan_idx );
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


void remove_reference_phase( cuDoubleComplex *J, cuDoubleComplex *Jref )
{
    cuDoubleComplex PP0norm, PQ0norm, QP0norm, QQ0norm;
    double PPscale = 1.0/cuCabs( Jref[0] ); // = 1/|PP|
    double PQscale = 1.0/cuCabs( Jref[1] ); // = 1/|PQ|
    double QPscale = 1.0/cuCabs( Jref[2] ); // = 1/|QP|
    double QQscale = 1.0/cuCabs( Jref[3] ); // = 1/|QQ|

    PP0norm = make_cuDoubleComplex( PPscale*cuCreal(Jref[0]), PPscale*cuCimag(Jref[0]) ); // = PP/|PP|
    PQ0norm = make_cuDoubleComplex( PQscale*cuCreal(Jref[1]), PQscale*cuCimag(Jref[1]) ); // = PQ/|PQ|
    QP0norm = make_cuDoubleComplex( QPscale*cuCreal(Jref[2]), QPscale*cuCimag(Jref[2]) ); // = QP/|QP|
    QQ0norm = make_cuDoubleComplex( QQscale*cuCreal(Jref[3]), QQscale*cuCimag(Jref[3]) ); // = QQ/|QQ|

    // Essentially phase rotations
    if (isfinite(PPscale)) { J[0] = cuCdiv( J[0], PP0norm ); }
    if (isfinite(PQscale)) { J[1] = cuCdiv( J[1], PQ0norm ); }
    if (isfinite(QPscale)) { J[2] = cuCdiv( J[2], QP0norm ); }
    if (isfinite(QQscale)) { J[3] = cuCdiv( J[3], QQ0norm ); }
}

void zero_PQ_and_QP( cuDoubleComplex *J )
/* For J = [ PP, PQ ], set PQ and QP to 0 for all antennas
 *         [ QP, QQ ]
 */
{
    J[1] = make_cuDoubleComplex( 0.0, 0.0 );
    J[2] = make_cuDoubleComplex( 0.0, 0.0 );
}

void parse_calibration_correction_file( uint32_t gpstime, struct calibration *cal )
/* Retrieve the PQ phase correction from pq_phase_correction.txt for the given
 * gps time
 */
{
    char buffer[CAL_BUFSIZE];
    sprintf( buffer, "%s/pq_phase_correction.txt", RUNTIME_DIR );
    FILE *f = fopen( buffer, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error:: Cannot find file '%s'\n", buffer );
        exit(EXIT_FAILURE);
    }

    // Values to be read in
    uint32_t from, to;
    double slope = 0.0, slope_tmp; // rad/Hz
    double offset = 0.0, offset_tmp; // rad
    char ref_tile_name[32];
    bool found = false;

    // Read in the file line by line
    while (fgets( buffer, CAL_BUFSIZE, f ) != NULL)
    {
        // Only parse if the first character is numeric
        if (!(buffer[0] >= '0' && buffer[0] <= '9'))
            continue;

        // Parse the line
        if (sscanf( buffer, "%u %u %lf %lf %s", &from, &to, &slope_tmp, &offset_tmp, ref_tile_name ) != 5)
        {
            fprintf( stderr, "error: parse_calibration_correction_file: cannot parse line \"%s\"",
                    buffer );
            exit(EXIT_FAILURE);
        }

        // If the given gpstime is in the right range, keep the values and
        // stop looking
        if (from <= gpstime && gpstime <= to)
        {
            slope = slope_tmp;
            offset = offset_tmp;
            found = true;
            break;
        }
    }

    // Close the file
    fclose( f );

    // Install the read-in values into the calibration struct, as appropriate
    if (cal->ref_ant == NULL) // i.e. the user hasn't opted to override using a custom reference tile
    {
        if (found)
        {
            cal->ref_ant = (char *)malloc( strlen(ref_tile_name) + 1 );
            strcpy( cal->ref_ant, ref_tile_name );
        }
        else
        {
            cal->ref_ant = (char *)malloc( strlen("NONE") + 1 );
            strcpy( cal->ref_ant, "NONE" );
        }
    }

    if (cal->custom_pq_correction == false) // i.e. the user has not provided a custom phase correction
    {
        if (found)
        {
            cal->phase_slope  = slope;
            cal->phase_offset = offset;
        }
    }
}

void apply_calibration_corrections( struct calibration *cal, cuDoubleComplex *D, MetafitsMetadata *obs_metadata,
        int coarse_chan_idx, logger *log )
/* Optionally apply certain corrections/adjustments to the calibration solutions
 * given in D. The corrections to be applied are stored in the calibration struct CAL.
 * Three corrections are applied here:
 *   (1) Subtract the phases from each antenna by a reference antenna
 *   (2) Remove the PQ and QP (i.e. off-diagonal) terms
 *   (3) Apply a phase slope to the QQ terms
 */
{
    char log_message[256];

    // Three locally defined booleans for whether to do the corrections
    bool apply_ref_ant         = true; // Whether this should be true is checked below
    bool apply_zero_PQ_and_QP  = !(cal->keep_cross_terms);
    bool apply_phase_slope     = (cal->phase_offset != 0.0 || cal->phase_slope != 0.0);

    Antenna *Aref = NULL; // The reference antenna struct
    if (strcmp(cal->ref_ant, "NONE") == 0) // i.e. the user explicitly wanted to NOT apply the reference antenna correction
    {
        apply_ref_ant = false;
    }
    else
    {
        Aref = find_antenna_by_name( obs_metadata, cal->ref_ant );

        // If the lookup failed, then quit with an error
        if (Aref == NULL)
        {
            fprintf( stderr, "error: apply_calibration_corrections: "
                    "tile %s does not exist in this observation\n", cal->ref_ant );
            exit(EXIT_FAILURE);
        }
    }

    // Quick check: if all three corrections are turned off, then abort here
    if (!apply_ref_ant && !apply_zero_PQ_and_QP && !apply_phase_slope)
        return;

    // Report on what's being done, if requested
    if (log)
    {
        if (apply_ref_ant)
        {
            sprintf( log_message, "Using tile %s as reference tile for rotating phases", cal->ref_ant );
            logger_timed_message( log, log_message );
        }
        else
        {
            logger_timed_message( log, "Rotation of phases relative to reference tile is turned OFF" );
        }

        if (apply_zero_PQ_and_QP)
            logger_timed_message( log, "Setting off-diagonal terms of calibration Jones matrix (PQ and QP) to zero" );
        else
            logger_timed_message( log, "Retaining off-diagonal terms of calibration Jones matrix (PQ and QP)" );

        if (apply_phase_slope)
        {
            sprintf( log_message, "Applying phase slope %e*FREQ + %e (rad) to PP", cal->phase_slope, cal->phase_offset );
            logger_timed_message( log, log_message );
        }
        else
            logger_timed_message( log, "Phase slope correction turned OFF" );
    }

    // Variables for converting slope and offset to a complex phase
    //     z = exp(i*phi) = cos(phi) + i*sin(phi),
    // where
    //     phi = slope*freq + offset

    double          phi; // rad
    cuDoubleComplex z;   // complex phase

    long int freq_ch; // Hz
    long int frequency  = obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_start_hz; // Hz
    int      chan_width = obs_metadata->corr_fine_chan_width_hz; // Hz

    uintptr_t d_idx, dref_idx;

    uintptr_t nant    = obs_metadata->num_ants;
    uintptr_t nantpol = obs_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t nchan   = obs_metadata->num_volt_fine_chans_per_coarse;

    // A temporary copy of the reference antenna matrix, so that we don't clobber it midway
    // through this operation!
    cuDoubleComplex Dref[nantpol*nantpol];

    // Go through all the antennas and divide the phases by the reference antenna
    uintptr_t ant, ch;
    for (ch = 0; ch < nchan; ch++)
    {
        if (apply_ref_ant)
        {
            // Make a copy of the reference Jones matrix for this channel
            // reference antenna and channel
            dref_idx = J_IDX(Aref->ant,ch,0,0,nchan,nantpol);
            cp2x2( &(D[dref_idx]), Dref );
        }

        if (apply_phase_slope)
        {
            // Convert the slope and offset into a complex phase
            freq_ch = frequency + ch*chan_width;                // The frequency of this fine channel (Hz)
            phi = cal->phase_slope*freq_ch + cal->phase_offset; // (rad)
            z = make_cuDoubleComplex( cos(phi), sin(phi) );
        }

        for (ant = 0; ant < nant; ant++)
        {
            // A pointer to the (first element of) the Jones matrix for this
            // antenna and channel
            d_idx = J_IDX(ant,ch,0,0,nchan,nantpol);

            // Divide through a reference antenna...
            if (apply_ref_ant)
            {
                remove_reference_phase( &(D[d_idx]), Dref );
//if (ch == 0) fprintf( stderr, "Dividing antenna %lu\n", ant );
            }

            // ...zero the off-diagonal terms...
            if (apply_zero_PQ_and_QP)
                zero_PQ_and_QP( &(D[d_idx]) );

            // ...and apply the phase correction to the PP element (pol 0,0)
            if (apply_phase_slope)
            {
                d_idx = J_IDX(ant,ch,0,0,nchan,nantpol);
                D[d_idx] = cuCmul( D[d_idx], z );
            }
        }
    }
}

void parse_flagged_tilenames_file( char *filename, struct calibration *cal )
{
    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error: parse_flagged_tilenames_file: "
                "could not open file '%s' for reading. Exiting\n", filename );
        exit(EXIT_FAILURE);
    }

    // Count the number of tiles in the list
    char tilename[16]; // More than enough space for tilenames
    cal->nflags = 0;
    while (fscanf( f, "%s", tilename ) != EOF)
        cal->nflags++;

    // Allocate memory for the list of tilenames in the calibration struct
    cal->flagged_tilenames = (char **)malloc( cal->nflags * sizeof(char *) );
    int i;
    for (i = 0; i < cal->nflags; i++)
        cal->flagged_tilenames[i] = (char *)malloc( 16 * sizeof(char) );

    // Rewind to the beginning of the file, and this time read them into
    // the calibration struct
    rewind( f );
    for (i = 0; i < cal->nflags; i++)
        fscanf( f, "%s", cal->flagged_tilenames[i] );

    // Close the file
    fclose( f );
}

bool tilename_is_flagged( char *tilename, struct calibration *cal )
{
    int i;
    for (i = 0; i < cal->nflags; i++)
    {
        if (strcmp( tilename, cal->flagged_tilenames[i] ) == 0) // A match!
        {
            return true;
        }
    }

    return false;
}

void set_flagged_tiles_to_zero( struct calibration *cal, MetafitsMetadata *obs_metadata, cuDoubleComplex *D )
{
    // Create a "zero" matrix that will be copied
    cuDoubleComplex Zero[4];
    int i;
    for (i = 0; i < 4; i++)
        Zero[i] = make_cuDoubleComplex( 0.0, 0.0 );

    // Loop through the tilenames listed in the calibration struct
    char     *tilename; // The tilename in question
    Antenna  *Ant;      // The Antenna struct for a given tilename
    uint32_t  ant;      // The corresponding antenna number
    uintptr_t nchan   = obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t nantpol = obs_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t ch;       // A particular fine channel
    uintptr_t d_idx;    // Idx into the D array

    for (i = 0; i < cal->nflags; i++)
    {
        // Pointer to the tilename in question
        tilename = cal->flagged_tilenames[i];

        if (!tilename_is_flagged( tilename, cal ))
            continue;

        // Find the corresponding antenna in the observation
        Ant = find_antenna_by_name( obs_metadata, tilename );
        ant = Ant->ant;

        // Loop through the fine channels
        for (ch = 0; ch < nchan; ch++)
        {
            // Get the index into the D array
            d_idx = J_IDX(ant,ch,0,0,nchan,nantpol);

            // Set it to zero
            cp2x2( Zero, &(D[d_idx]) );
        }
    }
}

void init_calibration( struct calibration *cal )
{
    cal->caldir            = NULL;
    cal->ref_ant           = NULL;
    cal->nflags            = 0;
    cal->flagged_tilenames = NULL;
}

void free_calibration( struct calibration *cal )
{
    if (cal->caldir != NULL)
        free( cal->caldir );

    if (cal->ref_ant != NULL)
        free( cal->ref_ant );

    int i;
    for (i = 0; i < cal->nflags; i++)
    {
        free( cal->flagged_tilenames[i] );
    }
    free( cal->flagged_tilenames );
}
