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

/**
 * vmLoadRTSSolution
 * =================
 *
 * Read in the RTS solution from the DI_Jones and Bandpass files in the
 * VM->CAL.CALDIR directory. It only reads the solution for the coarse channel
 * specified in VM->COARSE_CHAN_IDX. The output is a set of Jones matrices (D)
 * for each antenna and (non-flagged) fine channel.
 *
 * The Jones matrices (D) are output in the (Q,P) basis, i.e.
 * [ 1  2 ] = [ QQ  QP ]
 * [ 3  4 ] = [ PQ  PP ]
 *
 * The RTS bandpass solutions are used if VM->CAL.USE_BANDPASS is set,
 * otherwise only the coarse channel solutions are used. (This has no effect
 * for Offringa-style solutions.)
 *
 * This function assumes that a buffer for the D matrices has already been
 * allocated.
 */
void vmLoadRTSSolution( vcsbeam_context *vm )
{
    // Shorthand variables
    bool        use_bandpass = vm->cal.use_bandpass;
    const char *caldir       = vm->cal.caldir;

    // Find the "GPUBox" number for this coarse channel
    uintptr_t gpubox_number = vm->cal_metadata->metafits_coarse_chans[vm->coarse_chan_idx].corr_chan_number + 1;

    // With the gpubox number in hand, construct the filenames for the
    // DI_Jones and Bandpass files
    char dijones_path[CAL_BUFSIZE];
    char bandpass_path[CAL_BUFSIZE];

    sprintf( dijones_path,  "%s/DI_JonesMatrices_node%03lu.dat", caldir, gpubox_number );
    sprintf( bandpass_path, "%s/BandpassCalibration_node%03lu.dat", caldir, gpubox_number );

    // Allocate memory for the Jones arrays
    uintptr_t nant    = vm->cal_metadata->num_ants;
    uintptr_t ninputs = vm->cal_metadata->num_rf_inputs;
    uintptr_t nvispol = vm->cal_metadata->num_visibility_pols; // = 4 (PP, PQ, QP, QQ)
    uintptr_t nantpol = vm->cal_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t nchan   = vm->cal_metadata->num_corr_fine_chans_per_coarse;
    uintptr_t vcs_nchan = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t interp_factor = vcs_nchan / nchan;
    uintptr_t ant, ch; // For loop variables

    // Write out the channel --> gpubox number mapping
    if (vm->log)
    {
        sprintf( vm->log_message, "Receiver channel #%lu --> GPUBox #%lu",
                vm->cal_metadata->metafits_coarse_chans[vm->coarse_chan_idx].rec_chan_number,
                gpubox_number );
        logger_timed_message( vm->log, vm->log_message );
    }

    // Temporary arrays for DI Jones matrices ('Dd') and Bandpass matrices ('Db')
    cuDoubleComplex  **Dd = (cuDoubleComplex ** )calloc( nant, sizeof(cuDoubleComplex * ) );
    cuDoubleComplex ***Db = (cuDoubleComplex ***)calloc( nant, sizeof(cuDoubleComplex **) );

    // The array for the final output product (D = Dd x Db)
    // This will have the same dimensions as the final Jones matrix, so use
    // J_IDX for indexing
    cuDoubleComplex A[nvispol];

    for (ant = 0; ant < nant; ant++)
    {
        Dd[ant] = (cuDoubleComplex * )calloc( nvispol, sizeof(cuDoubleComplex  ) );
        Db[ant] = (cuDoubleComplex **)calloc( nchan,   sizeof(cuDoubleComplex *) );

        for (ch = 0; ch < nchan; ch++)
            Db[ant][ch] = (cuDoubleComplex *)calloc( nvispol, sizeof(cuDoubleComplex) );
    }

    // Read in the DI Jones file
    read_dijones_file(Dd, A, NULL, vm->cal_metadata->num_ants, dijones_path);

    // Invert the alignment matrix and multiply it to the gains matrix
    // Eq. (29) in Ord et al. (2019)
    cuDoubleComplex Ainv[nvispol];
    inv2x2( A, Ainv );

    // Read in the Bandpass file
    if (use_bandpass)
        read_bandpass_file( NULL, Db, vm->cal_metadata, bandpass_path );

    // Make the master mpi thread print out the antenna names of both
    // obs and cal metafits. "Header" printed here, actual numbers
    // printed inside antenna for loop below
    if (vm->log)
    {
        if (vm->log->world_rank == 0)
        {
            logger_message( vm->log, "" );
            logger_timed_message( vm->log, "Calibration metafits info:" );
            logger_timed_message( vm->log, "--------------------------" );
            logger_timed_message( vm->log, "|Input| Ant |Tile ID| TileName |Pol|VCS order|" );

            uintptr_t i;
            for (i = 0; i < ninputs; i++)
            {
                sprintf( vm->log_message, "| %3u | %3u | %5u | %8s | %c |  %3u   |",
                        vm->cal_metadata->rf_inputs[i].input,
                        vm->cal_metadata->rf_inputs[i].ant,
                        vm->cal_metadata->rf_inputs[i].tile_id,
                        vm->cal_metadata->rf_inputs[i].tile_name,
                        *(vm->cal_metadata->rf_inputs[i].pol),
                        vm->cal_metadata->rf_inputs[i].vcs_order
                       );
                logger_timed_message( vm->log, vm->log_message );
            }

            logger_message( vm->log, "" );
        }
    }

    // Form the "fine channel" DI gain (the "D" in Eqs. (28-30), Ord et al. (2019))
    // Do this for each of the _voltage_ observation's fine channels (use
    // nearest-neighbour interpolation). This is "applying the bandpass" corrections.
    uintptr_t cal_ch, obs_ant, dd_idx;
    uintptr_t d_idx;
    uintptr_t i; // (i)dx into rf_inputs
    Rfinput *obs_rfinput, *cal_rfinput;
    for (i = 0; i < ninputs; i++)
    {
        // Order the Jones matrices by the TARGET OBSERVATION metadata
        // (not the calibration metadata)
        obs_rfinput = &(vm->obs_metadata->rf_inputs[i]);

        // Only have to loop once per tile, so skip the 'Y's
        if (*(obs_rfinput->pol) == 'Y')
            continue;

        // Match up antenna indices between the cal and obs metadata
        cal_rfinput = find_matching_rf_input( vm->cal_metadata, obs_rfinput );

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
            if (use_bandpass)
                mult2x2d( Dd[dd_idx], Db[dd_idx][cal_ch], &(vm->D[d_idx]) );
            else
                cp2x2( Dd[dd_idx], &(vm->D[d_idx]) );

            // Multiply in (the inverse of) the alignment matrix
            mult2x2d( &(vm->D[d_idx]), Ainv, &(vm->D[d_idx]) );

            // The RTS matrices are apparently in (p,q)<-(p,q) basis. The
            // following converts to (q,p)<-(q,p)
            reverse2x2( &(vm->D[d_idx]), &(vm->D[d_idx]) );
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
}

void read_dijones_file( cuDoubleComplex **Dd, cuDoubleComplex *A, double *amp, uintptr_t nant, char *fname )
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
    double *Adbl = (double *)A; // The alignment matrix, cast as a double array for easy reading
    double J[NDBL_PER_JONES]; // "Raw" Jones matrix

    // Read in the amplitude (SM: I have no idea what this value represents)
    fscanf( fp, "%lf", &val );
    if (amp != NULL)
        *amp = val;

    // Read in the alignment ("A")
    int i;
    for (i = 0; i < NDBL_PER_JONES; i++)
        fscanf( fp, "%lf,", &Adbl[i] );

    // Finally, read in the Jones ("J") matrices
    uintptr_t ant;
    for (ant = 0; ant < nant; ant++)
    {
        for (i = 0; i < NDBL_PER_JONES; i++)
            fscanf( fp, "%lf,", &J[i] );

        cp2x2( (cuDoubleComplex *)J, Dd[ant] );
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


/**
 * vmLoadOffringaSolution
 * ======================
 *
 * Offringa's own documentation for the format of these binary files (copied from an email
 * dated 4 May 2016 from franz.kirsten@curtin.edu.au to sammy.mcsweeney@gmail.com). This is
 * itself copied from Andre Offringa's original C++ code, in his Anoko repository, in
 * mwa-reduce/solutionfile.h:
 *
 * The solution file is used for storing calibration Jones matrices.
 * The format is as follows:
  *  Bytes |  Description
  * -------+---------------
  *  0- 7  |  string intro ; 8-byte null terminated string "MWAOCAL"
  *  8-11  |  int fileType ; always 0, reserved for indicating something
  *           other than complex Jones solutions
  * 12-15  |  int structureType ; always 0, reserved for indicating
  *           different ordering
  * 16-19  |  int intervalCount ; Number of solution intervals in file
  * 20-23  |  int antennaCount ; Number of antennas that were in the
  *           measurement set (but were not necessary all solved for)
  * 24-27  |  int channelCount ; Number of channels in the measurement set
  * 28-31  |  int polarizationCount ; Number of polarizations solved for
  *           -- always four.
  * 32-39  |  double startTime ; Start time of solutions (AIPS time)
  * 40-47  |  double endTime ; End time of solutions (AIPS time)
  * -------+-------------------
  * After the header follow 2 x nSolution doubles, with
  *
  * nSolutions = nIntervals * nAntennas * nChannels * nPols
  *
  * Ordered in the way as given, so:
  * double 0 : real of interval 0, antenna 0, channel 0, pol 0
  * double 1 : imaginary of interval 0, antenna 0, channel 0, pol 0
  * double 2 : real of interval 0, antenna 0, channel 0, pol 1
  * ...
  * double 8 : real of interval 0, antenna 0, channel 1, pol 0
  * double nChannel x 8 : real of interval 0, antenna 1, channel 0, pol 0
  * etc.
  *
  * ints are here always 32 bits unsigned integers, doubles are IEEE
  * double precision 64 bit floating points.
  * If a solution is not available, either because its data no data were
  * selected during calibration for this interval
  * or because the calibration diverged, a "NaN" will be stored in the
  * doubles belonging to that solution.
  *
  * ----------------------------------------------
  *
  * The input file is taken from VM->CAL.CALDIR.
  *
  * Only the coarse channel specified in VM->COARSE_CHAN_IDX will be read in.
  *
  * This function assumes that memory for the Jones matrices (D) has already
  * been allocated.
  */
void vmLoadOffringaSolution( vcsbeam_context *vm )
{
    // Shorthand variables
    int coarse_chan_idx = vm->coarse_chan_idx;

    // Open the calibration file for reading
    FILE *fp = NULL;
    fp = fopen(vm->cal.caldir,"r");
    if (fp == NULL)
    {
        fprintf( stderr, "Failed to open %s: quitting\n", vm->cal.caldir );
        exit(EXIT_FAILURE);
    }

    // Read in the necessary information from the header

    uint32_t intervalCount, antennaCount, channelCount, polarizationCount;
    uint32_t nant   = vm->cal_metadata->num_ants;
    uint32_t nchan  = vm->cal_metadata->num_corr_fine_chans_per_coarse;
    uint32_t nChan  = nchan * vm->cal_metadata->num_metafits_coarse_chans;
    uint32_t ninput = vm->cal_metadata->num_rf_inputs;
    uintptr_t nantpol = vm->cal_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t vcs_nchan = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t interp_factor = vcs_nchan / nchan;
    uintptr_t nvispol = vm->cal_metadata->num_visibility_pols; // = 4 (PP, PQ, QP, QQ)

    // Make another dummy matrix for reading in
    cuDoubleComplex Dread[nvispol];

    fseek(fp, 16, SEEK_SET);
    fread(&intervalCount,     sizeof(uint32_t), 1, fp);
    fread(&antennaCount,      sizeof(uint32_t), 1, fp);
    fread(&channelCount,      sizeof(uint32_t), 1, fp);
    fread(&polarizationCount, sizeof(uint32_t), 1, fp);

    // Error-checking the info extracted from the header,
    // making sure it matches the metadata
    if (intervalCount > 1)
    {
        fprintf( stderr, "Warning: Only the first interval in the calibration " );
        fprintf( stderr, "solution (%s) will be used\n", vm->cal.caldir );
    }
    if (antennaCount != nant)
    {
        fprintf( stderr, "Error: Calibration solution (%s) ", vm->cal.caldir );
        fprintf( stderr, "contains a different number of antennas (%d) ", antennaCount );
        fprintf( stderr, "than specified (%d)\n", nant );
        exit(EXIT_FAILURE);
    }
    if (channelCount != nChan)
    {
        fprintf( stderr, "Warning: Calibration solution (%s) ", vm->cal.caldir );
        fprintf( stderr, "contains a different number (%d) ", channelCount );
        fprintf( stderr, "than the expected (%d) channels.\n", nChan );
        nChan = channelCount;
        nchan = nChan / vm->cal_metadata->num_metafits_coarse_chans;
        interp_factor = vcs_nchan / nchan;
        fprintf( stderr, "Assuming calibration channels are "
                "%d kHz\n", vm->cal_metadata->coarse_chan_width_hz / nchan / 1000 );
    }
    if (coarse_chan_idx >= (int)vm->cal_metadata->num_metafits_coarse_chans)
    {
        fprintf( stderr, "Error: Requested coarse channel number (%d) ", coarse_chan_idx );
        fprintf( stderr, "is more than the number of channels " );
        fprintf( stderr, "available in the calibration solution (%s)\n", vm->cal.caldir );
        exit(EXIT_FAILURE);
    }

    // Iterate through antennas and channels
    uint32_t ant; // The antenna number (as defined in metafits "Antenna")
    uint32_t i;   // Just a dummy index into the metafits array of rf_inputs
    uint32_t ch;  // (Fine) channel number (within given coarse channel)
    uint32_t Ch;  // (Fine) channel number (within total set of channels)
    uint32_t vch; // (Fine) channel number (within coarse channel in voltage obs)
    uint32_t Interval = 0;  // Only use the first interval
    long fpos; // Position within the gains file
    uintptr_t d_idx; // Idx into the D array
    Rfinput *rfinput = NULL;

    for (i = 0; i < ninput; i++)
    {
        // We only need information for each antenna, so skip one of the pols
        rfinput = &(vm->cal_metadata->rf_inputs[i]);
        if (*(rfinput->pol) == 'Y')
            continue;

        // Get the antenna number
        ant = rfinput->ant;

        // Loop over channels
        for (ch = 0; ch < nchan; ch++)
        {
            // Translate from "fine channel number within coarse channel"
            // to "fine channel number within whole observation"
            Ch = ch + coarse_chan_idx*nchan;

            // Move the file pointer to the correct place
            fpos = OFFRINGA_HEADER_SIZE_BYTES +
                Interval * (nant * nChan * JONES_SIZE_BYTES) +
                ant      *        (nChan * JONES_SIZE_BYTES) +
                Ch       *                (JONES_SIZE_BYTES);
            fseek( fp, fpos, SEEK_SET );

            // Read in one Jones matrix
            fread( Dread, 1, JONES_SIZE_BYTES, fp );

            // If there are any nans, set them to zero instead
            // Assume that if there are any nans in the Jones matrix, then
            // EVERY element in the Jones matrix is a nan. Therefore, only need
            // to check one element.
            if (isnan(Dread[0].x))
                memset( Dread, 0, JONES_SIZE_BYTES );

            // Copy this matrix into every corresponding "voltage" fine channel
            for (vch = ch*interp_factor; vch < (ch + 1)*interp_factor; vch++)
            {
                // Get the destination index
                d_idx = J_IDX(ant,vch,0,0,vcs_nchan,nantpol);

                // Copy it across
                cp2x2( Dread, &(vm->D[d_idx]) );
            }
        }
    }

    // Close the file
    fclose(fp);
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
    if (isfinite(PPscale))  J[0] = cuCdiv( J[0], PP0norm );
    if (isfinite(PQscale))  J[1] = cuCdiv( J[1], PQ0norm );
    if (isfinite(QPscale))  J[2] = cuCdiv( J[2], QP0norm );
    if (isfinite(QQscale))  J[3] = cuCdiv( J[3], QQ0norm );
}

void zero_PQ_and_QP( cuDoubleComplex *J )
/* For J = [ QQ, QP ], set PQ and QP to 0 for all antennas
 *         [ PQ, PP ]
 */
{
    J[1] = make_cuDoubleComplex( 0.0, 0.0 );
    J[2] = make_cuDoubleComplex( 0.0, 0.0 );
}

void parse_calibration_correction_file( uint32_t gpstime, calibration *cal )
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

void apply_calibration_corrections( calibration *cal, cuDoubleComplex *D, MetafitsMetadata *obs_metadata,
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
            }

            // ...zero the off-diagonal terms...
            if (apply_zero_PQ_and_QP)
                zero_PQ_and_QP( &(D[d_idx]) );

            // ...and apply the phase correction:
            // DZ = [ d00 d01 ] [ z 0 ]
            //      [ d10 d11 ] [ 0 1 ]
            //
            //    = [ d00*z  d01 ]
            //      [ d10*z  d11 ]
            if (apply_phase_slope)
            {
                d_idx = J_IDX(ant,ch,0,0,nchan,nantpol);
                D[d_idx] = cuCmul( D[d_idx], z );

                d_idx = J_IDX(ant,ch,1,0,nchan,nantpol);
                D[d_idx] = cuCmul( D[d_idx], z );
            }
        }
    }
}

void vmParseFlaggedTilenamesFile( char *filename, calibration *cal )
{
    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error: vmParseFlaggedTilenamesFile: "
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

bool tilename_is_flagged( char *tilename, calibration *cal )
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

void vmSetCustomTileFlags( vcsbeam_context *vm, char *filename, calibration *cal )
{
    // If no filename given, just count the number of active tiles and return
    if (filename == NULL)
    {
        vmSetNumNotFlaggedRFInputs( vm );
        return;
    }

    // Parse the given file for tilenames
    vmParseFlaggedTilenamesFile( filename, cal );

    // Create a "zero" matrix that will be copied
    cuDoubleComplex Zero[4];
    int i;
    for (i = 0; i < 4; i++)
        Zero[i] = make_cuDoubleComplex( 0.0, 0.0 );

    // Loop through the tilenames listed in the calibration struct
    char     *tilename; // The tilename in question
    Antenna  *Ant;      // The Antenna struct for a given tilename
    uint32_t  ant;      // The corresponding antenna number
    uintptr_t nchan   = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t nantpol = vm->obs_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t ch;       // A particular fine channel
    uintptr_t d_idx;    // Idx into the D array

    for (i = 0; i < cal->nflags; i++)
    {
        // Pointer to the tilename in question
        tilename = cal->flagged_tilenames[i];

        if (!tilename_is_flagged( tilename, cal ))
            continue;

        // Find the corresponding antenna in the observation
        Ant = find_antenna_by_name( vm->obs_metadata, tilename );
        ant = Ant->ant;

        // Loop through the fine channels
        for (ch = 0; ch < nchan; ch++)
        {
            // Get the index into the D array
            d_idx = J_IDX(ant,ch,0,0,nchan,nantpol);

            // Set it to zero
            cp2x2( Zero, &(vm->D[d_idx]) );
        }
    }

    // Count the number of active tiles
    vmSetNumNotFlaggedRFInputs( vm );
}

void init_calibration( calibration *cal )
{
    cal->caldir            = NULL;
    cal->ref_ant           = NULL;
    cal->nflags            = 0;
    cal->flagged_tilenames = NULL;
}

void free_calibration( calibration *cal )
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
