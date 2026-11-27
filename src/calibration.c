/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <fitsio.h>
#include <mwalib.h>

#include "vcsbeam.h"
#include "gpu_macros.h"

/**
 * Reads in a calibration solution from an Offringa-style solution file.
 *
 * A full description of the Offringa file format can be found
 * [here](@ref offringafileformat).
 *
 * The input file is taken from `vm&rarr;cal.caldir`.
 *
 * Only the coarse channel specified in `vm&rarr;coarse_chan_idx` will be read
 * in.
 *
 * This function assumes that memory for the Jones matrices (D) has already
 * been allocated.
 */
void vmLoadOffringaSolution( vcsbeam_context *vm )
{
    // Shorthand variables
    int coarse_chan_idx = vm->cal_coarse_chan_idxs_to_process[0];

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
    uint32_t ninput = vm->cal_metadata->num_rf_inputs;
    uintptr_t nantpol = vm->cal_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t vcs_nchan = vm->nfine_chan;
    uintptr_t nvispol = vm->cal_metadata->num_visibility_pols; // = 4 (PP, PQ, QP, QQ)

    // Set these values later after checking the calibration solution size
    uint32_t nchan;
    uint32_t nChan;
    uintptr_t interp_factor;

    // Make another dummy matrix for reading in
    gpuDoubleComplex Dread[nvispol];

    fseek(fp, 16, SEEK_SET);
    fread(&intervalCount,     sizeof(uint32_t), 1, fp);
    fread(&antennaCount,      sizeof(uint32_t), 1, fp);
    fread(&channelCount,      sizeof(uint32_t), 1, fp);
    fread(&polarizationCount, sizeof(uint32_t), 1, fp);

#ifdef DEBUG
    fprintf( stderr, "*** Reading Offringa-style solution, DEBUG ***\n" );
    fprintf( stderr, "Quantities determined via vcsbeam context struct, DEBUG\n" );
    fprintf( stderr, "nant = vm->cal_metadata->num_ants = %u\n", nant );
    fprintf( stderr, "nchan = vm->cal_metadata->num_corr_fine_chans_per_coarse = %u\n", nchan );
    fprintf( stderr, "nChan = nchan * vm->mpi_size = %u\n", nChan );
    fprintf( stderr, "ninput = vm->cal_metadata->num_rf_inputs = %u\n", ninput );
    fprintf( stderr, "nantpol = vm->cal_metadata->num_ant_pols = %lu (should be 2)\n", nantpol );
    fprintf( stderr, "vcs_nchan = vm->nfine_chan = %lu\n", vcs_nchan );
    fprintf( stderr, "interp_factor = vcs_nchan / nchan = %lu\n", interp_factor );
    fprintf( stderr, "nvispol = vm->cal_metadata->num_visibility_pols = %lu (should be 4)\n", nvispol );

    fprintf( stderr, "Quantities read from binary file, DEBUG\n" );
    fprintf( stderr, "intervalCount = %u\n", intervalCount );
    fprintf( stderr, "antennaCount = %u\n", antennaCount );
    fprintf( stderr, "channelCount = %u\n", channelCount );
    fprintf( stderr, "polarizationCount = %u\n", polarizationCount );
#endif

    // Error-checking the info extracted from the header,
    // making sure it matches the metadata
    if (intervalCount > 1)
    {
        snprintf( vm->log_message, sizeof(vm->log_message),
                "Warning: Only the first interval in the calibration "
                "solution (%s) will be used", vm->cal.caldir );
        logger_timed_message( vm->log, vm->log_message );
    }
    if (antennaCount != nant)
    {
        fprintf( stderr, "Error: Calibration solution (%s) ", vm->cal.caldir );
        fprintf( stderr, "contains a different number of antennas (%u) ", antennaCount );
        fprintf( stderr, "than specified (%u)\n", nant );
        exit(EXIT_FAILURE);
    }

    // Now set the calibration solution settings correctly
    snprintf( vm->log_message, sizeof(vm->log_message),
            "Calibration solution has (%u) total channels", channelCount );
    logger_timed_message( vm->log, vm->log_message );
    snprintf( vm->log_message, sizeof(vm->log_message),
            "    Calibration solution file: %s", vm->cal.caldir );
    logger_timed_message( vm->log, vm->log_message );
    nChan = channelCount;

    // Figure out the correct number of fine channels in a coarse channel from 
    // the calibration solution
    if (vm->cal.picket_fence)
    {
        // Assuming calibration solution has same number of coarse channels as 
        // the number of MPI processes for picket fenced data
        snprintf( vm->log_message, sizeof(vm->log_message),
                "Assuming %u coarse channels in the calibration "
                "solution for picket fence data", vm->mpi_size );
        logger_timed_message( vm->log, vm->log_message );
        nchan = nChan / vm->mpi_size; 
    }
    else
    {
        // Assuming 24 coarse channels in calibration solution for contiguous 
        // data
        snprintf( vm->log_message, sizeof(vm->log_message),
                "Assuming %d coarse channels in the calibration "
                "solution for contiguous data", 24 );
        logger_timed_message( vm->log, vm->log_message );
        nchan = nChan / 24;
    }
    interp_factor = vcs_nchan / nchan;
    snprintf( vm->log_message, sizeof(vm->log_message),
            "Assuming calibration channels are %d kHz",
            vm->cal_metadata->coarse_chan_width_hz / nchan / 1000 );
    logger_timed_message( vm->log, vm->log_message );
#ifdef DEBUG
    fprintf( stderr, "New nChan = %u\n", nChan );
    fprintf( stderr, "New nchan = %u\n", nchan );
    fprintf( stderr, "New interp_factor = %lu\n", interp_factor );
#endif
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
            Ch = ch + coarse_chan_idx * nchan;

            // Move the file pointer to the correct place
            fpos = OFFRINGA_HEADER_SIZE_BYTES +
                Interval * (nant * nChan * JONES_SIZE_BYTES) +
                ant      *        (nChan * JONES_SIZE_BYTES) +
                Ch       *                (JONES_SIZE_BYTES);
            fseek( fp, fpos, SEEK_SET );

//fprintf(stderr, "ftell=%lu    ", ftell(fp));
            // Read in one Jones matrix
            fread( Dread, 1, JONES_SIZE_BYTES, fp );

//fprintf(stderr, "Dread (Interval = %d, ant = %d, Ch = %d): ", Interval, ant, Ch); fprintf_complex_matrix( stderr, Dread );

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
                d_idx = D_IDX(ant,vch,0,0,vcs_nchan,nantpol);

                // Copy it across
                cp2x2( Dread, &(vm->D[d_idx]) );
            }
        }
    }

    // Close the file
    fclose(fp);
}

/**
 * Exits with a helpful message if a cfitsio call failed.
 */
static void check_fits_status( int status, const char *filename, const char *what )
{
    if (status == 0)
        return;
    char errtxt[FLEN_ERRMSG];
    fits_get_errstatus( status, errtxt );
    fprintf( stderr, "Error: %s in hyperdrive solution (%s): %s\n", what, filename, errtxt );
    exit(EXIT_FAILURE);
}

/**
 * Opens the first hyperdrive solution file in a list that covers the given
 * coarse channel.
 *
 * @param[in]  nfiles    The number of files in the list
 * @param[in]  files     The hyperdrive solution files
 * @param[in]  start_hz  Lower edge of the coarse channel (Hz)
 * @param[in]  end_hz    Upper edge of the coarse channel (Hz)
 * @param[out] fptr      The opened file (positioned at an arbitrary HDU)
 * @param[out] filename  The name of the file that was opened (to be freed by the caller)
 * @param[out] first     The index of the first solution channel within the coarse channel
 * @param[out] flags     Per-channel flags (nchan entries, to be freed by the caller)
 * @param[out] nchan     The number of solution channels within the coarse channel
 *
 * Each hyperdrive file covers a contiguous subband. The CHANBLOCKS table gives
 * the centre frequency of each solution channel, so the channels belonging to
 * the coarse channel are found by frequency. This makes no assumption about
 * how many coarse channels each file contains or where they start.
 */
static void open_hyperdrive_file_for_coarse_chan( int nfiles, char **files,
        double start_hz, double end_hz,
        fitsfile **fptr, char **filename, long *first, unsigned char **flags, long *nchan )
{
    for (int fi = 0; fi < nfiles; fi++)
    {
        char *f = files[fi];
        int status = 0;
        fitsfile *fp = NULL;
        fits_open_file( &fp, f, READONLY, &status );
        check_fits_status( status, f, "failed to open file" );

        fits_movnam_hdu( fp, BINARY_TBL, "CHANBLOCKS", 0, &status );
        check_fits_status( status, f, "no CHANBLOCKS HDU" );
        int col_freq = 0, col_flag = 0;
        long nrows = 0;
        fits_get_colnum( fp, CASEINSEN, "Freq", &col_freq, &status );
        fits_get_colnum( fp, CASEINSEN, "Flag", &col_flag, &status );
        fits_get_num_rows( fp, &nrows, &status );
        check_fits_status( status, f, "failed to read CHANBLOCKS table" );

        double *freq = (double *)malloc( nrows * sizeof(double) );
        unsigned char *flag = (unsigned char *)malloc( nrows );
        int anynul = 0;
        fits_read_col( fp, TDOUBLE, col_freq, 1, 1, nrows, NULL, freq, &anynul, &status );
        fits_read_col( fp, TBYTE,   col_flag, 1, 1, nrows, NULL, flag, &anynul, &status );
        check_fits_status( status, f, "failed to read CHANBLOCKS table" );

        long i0 = -1, n = 0;
        for (long i = 0; i < nrows; i++)
        {
            if (freq[i] >= start_hz && freq[i] < end_hz)
            {
                if (i0 < 0) i0 = i;
                n++;
            }
        }
        free( freq );

        if (n > 0)
        {
            *fptr = fp;
            *filename = strdup( f );
            *first = i0;
            *nchan = n;
            *flags = (unsigned char *)malloc( n );
            memcpy( *flags, flag + i0, n );
            free( flag );
            return;
        }

        free( flag );
        fits_close_file( fp, &status );
    }

    fprintf( stderr, "Error: None of the %d calibration solution file(s) cover the "
            "coarse channel %.3f - %.3f MHz\n", nfiles, start_hz/1e6, end_hz/1e6 );
    exit(EXIT_FAILURE);
}

/**
 * Reads in a calibration solution from native mwa_hyperdrive FITS file(s).
 *
 * The files are `vm->cal.cal_files` (see vmBindCalibrationFiles()), or just
 * `vm->cal.caldir` if no list was bound. Each hyperdrive file covers one 
 * contiguous subband (anything from one to all 24 coarse channels). 
 * The coarse channel to be processed by this MPI task
 * (`vm->cal_coarse_chan_idxs_to_process[0]`) is located by frequency, using the
 * `CHANBLOCKS` table, in the first file that contains it.
 *
 * Within a file, the relevant HDUs are:
 * - `SOLUTIONS`: double-precision image with FITS axes
 *   8 x nchan x ntiles x ntimeblocks. The 8 values are the real and imaginary
 *   parts of the Jones matrix elements in the order PP, PQ, QP, QQ.
 * - `TILES`: binary table with `Antenna` and `Flag` columns.
 * - `CHANBLOCKS`: binary table with `Freq` (Hz) and `Flag` columns.
 *
 * Only the first timeblock is used. NaN solutions, flagged tiles and flagged
 * channels are set to zero (i.e. flagged), as for vmLoadOffringaSolution().
 *
 * This function assumes that memory for the Jones matrices (D) has already
 * been allocated.
 */
void vmLoadHyperdriveSolution( vcsbeam_context *vm )
{
    int coarse_chan_idx = vm->cal_coarse_chan_idxs_to_process[0];
    uint32_t nant   = vm->cal_metadata->num_ants;
    uint32_t ninput = vm->cal_metadata->num_rf_inputs;
    uintptr_t nantpol   = vm->cal_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t vcs_nchan = vm->nfine_chan;
    uintptr_t nvispol   = vm->cal_metadata->num_visibility_pols; // = 4

    if (coarse_chan_idx >= (int)vm->cal_metadata->num_metafits_coarse_chans)
    {
        fprintf( stderr, "Error: Requested coarse channel number (%d) is more than "
                "the number of coarse channels in the calibration observation\n", coarse_chan_idx );
        exit(EXIT_FAILURE);
    }
    double chan_start_hz = vm->cal_metadata->metafits_coarse_chans[coarse_chan_idx].chan_start_hz;
    double chan_end_hz   = vm->cal_metadata->metafits_coarse_chans[coarse_chan_idx].chan_end_hz;

    // Find the file (and the range of channels in it) for this coarse channel
    int status = 0;
    fitsfile *fptr = NULL;
    char *filename = NULL;
    long first_chan, nchan;
    unsigned char *chan_flag = NULL;
    int nfiles = (vm->cal.ncal_files > 0) ? vm->cal.ncal_files : 1;
    char **files = (vm->cal.ncal_files > 0) ? vm->cal.cal_files : &(vm->cal.caldir);
    open_hyperdrive_file_for_coarse_chan( nfiles, files, chan_start_hz, chan_end_hz,
            &fptr, &filename, &first_chan, &chan_flag, &nchan );

    if (nchan > (long)vcs_nchan)
    {
        fprintf( stderr, "Error: Calibration solution (%s) has more channels per coarse "
                "channel (%ld) than the voltage data (%lu)\n", filename, nchan, vcs_nchan );
        exit(EXIT_FAILURE);
    }
    uintptr_t interp_factor = vcs_nchan / nchan;
    snprintf( vm->log_message, sizeof(vm->log_message),
            "Calibration solution: using %ld channels (%d kHz each) for "
            "coarse channel %d/%d", nchan,
            (int)(vm->cal_metadata->coarse_chan_width_hz / nchan / 1000),
            coarse_chan_idx, vm->ncoarse_chans );
    logger_timed_message( vm->log, vm->log_message );
    snprintf( vm->log_message, sizeof(vm->log_message),
            "    Calibration solution file: %s", filename );
    logger_timed_message( vm->log, vm->log_message );

    // --- Solutions image dimensions ---
    fits_movnam_hdu( fptr, IMAGE_HDU, "SOLUTIONS", 0, &status );
    check_fits_status( status, filename, "no SOLUTIONS HDU" );
    int naxis = 0;
    long naxes[4] = {0,0,0,0};
    fits_get_img_dim( fptr, &naxis, &status );
    check_fits_status( status, filename, "failed to read SOLUTIONS dimensions" );
    if (naxis != 4)
    {
        fprintf( stderr, "Error: SOLUTIONS HDU in %s has %d axes (expected 4)\n", filename, naxis );
        exit(EXIT_FAILURE);
    }
    fits_get_img_size( fptr, 4, naxes, &status );
    check_fits_status( status, filename, "failed to read SOLUTIONS dimensions" );

    long nChan = naxes[1], ntiles = naxes[2], ntimeblocks = naxes[3];
    if (naxes[0] != NDBL_PER_JONES)
    {
        fprintf( stderr, "Error: SOLUTIONS HDU in %s has %ld values per Jones matrix "
                "(expected %d)\n", filename, naxes[0], NDBL_PER_JONES );
        exit(EXIT_FAILURE);
    }
    if (ntiles != (long)nant)
    {
        fprintf( stderr, "Error: Calibration solution (%s) contains a different "
                "number of antennas (%ld) than specified (%u)\n", filename, ntiles, nant );
        exit(EXIT_FAILURE);
    }
    if (ntimeblocks > 1)
    {
        snprintf( vm->log_message, sizeof(vm->log_message),
                "Warning: Only the first timeblock in the calibration "
                "solution (%s) will be used", filename );
        logger_timed_message( vm->log, vm->log_message );
    }

    // --- Tiles table: map antenna number -> row in the solutions, and flags ---
    long *tile_row = (long *)malloc( nant * sizeof(long) );
    for (uint32_t a = 0; a < nant; a++)
        tile_row[a] = -1;
    short *tile_flag = (short *)calloc( ntiles, sizeof(short) );

    fits_movnam_hdu( fptr, BINARY_TBL, "TILES", 0, &status );
    check_fits_status( status, filename, "no TILES HDU" );
    int col_ant = 0, col_flag = 0;
    fits_get_colnum( fptr, CASEINSEN, "Antenna", &col_ant, &status );
    fits_get_colnum( fptr, CASEINSEN, "Flag", &col_flag, &status );
    check_fits_status( status, filename, "TILES table lacks Antenna/Flag columns" );
    for (long r = 0; r < ntiles; r++)
    {
        int antnum = 0, anynul = 0;
        fits_read_col( fptr, TINT,   col_ant,  r+1, 1, 1, NULL, &antnum,         &anynul, &status );
        fits_read_col( fptr, TSHORT, col_flag, r+1, 1, 1, NULL, &(tile_flag[r]), &anynul, &status );
        check_fits_status( status, filename, "failed to read TILES table" );
        if (antnum >= 0 && antnum < (int)nant)
            tile_row[antnum] = r;
    }

    // --- Read the Jones matrices ---
    fits_movnam_hdu( fptr, IMAGE_HDU, "SOLUTIONS", 0, &status );
    check_fits_status( status, filename, "no SOLUTIONS HDU" );

    double *buf = (double *)malloc( nchan * NDBL_PER_JONES * sizeof(double) );
    gpuDoubleComplex Dread[nvispol];
    const long timeblock = 0;

    for (uint32_t i = 0; i < ninput; i++)
    {
        Rfinput *rfinput = &(vm->cal_metadata->rf_inputs[i]);
        if (*(rfinput->pol) == 'Y')
            continue;

        uint32_t ant = rfinput->ant;
        if (ant >= nant || tile_row[ant] < 0)
        {
            fprintf( stderr, "Error: Antenna %u (%s) not found in the TILES table of %s\n",
                    ant, rfinput->tile_name, filename );
            exit(EXIT_FAILURE);
        }
        long row = tile_row[ant];

        // C-order index [timeblock][tile][chan][8], 1-indexed for cfitsio
        LONGLONG firstelem = (((LONGLONG)timeblock*ntiles + row)*nChan + first_chan)
                * NDBL_PER_JONES + 1;
        int anynul = 0;
        fits_read_img( fptr, TDOUBLE, firstelem, nchan*NDBL_PER_JONES, NULL, buf, &anynul, &status );
        check_fits_status( status, filename, "failed to read solutions" );

        for (long ch = 0; ch < nchan; ch++)
        {
            double *J = buf + ch*NDBL_PER_JONES;
            int bad = (tile_flag[row] != 0) || (chan_flag[ch] != 0);
            for (int k = 0; k < NDBL_PER_JONES; k++)
                if (isnan(J[k]))
                    bad = 1;

            for (uintptr_t p = 0; p < nvispol; p++)
                Dread[p] = bad ? make_gpuDoubleComplex( 0.0, 0.0 )
                               : make_gpuDoubleComplex( J[2*p], J[2*p+1] );

            for (uintptr_t vch = ch*interp_factor; vch < (uintptr_t)(ch + 1)*interp_factor; vch++)
            {
                uintptr_t d_idx = D_IDX(ant,vch,0,0,vcs_nchan,nantpol);
                cp2x2( Dread, &(vm->D[d_idx]) );
            }
        }
    }

    free( buf );
    free( tile_row );
    free( tile_flag );
    free( chan_flag );
    free( filename );
    fits_close_file( fptr, &status );
}

/**
 * Rotates the phases of the elements of a complex matrix with respect to a
 * reference matrix.
 *
 * @param[in,out] J     The matrix whose phases are to be rotated, \f${\bf J}\f$
 * @param[in]     Jref  The reference matrix, \f${\bf J}_\text{ref}\f$
 *
 * Each element in \f${\bf J}\f$ is divided by the corresponding element in
 * \f${\bf J}_\text{ref}\f$, normalised to unit magnitude.
 * That is, if
 * \f[
 * \begin{aligned}
 *     {\bf J} &= \begin{bmatrix} J_0 & J_1 \\ J_2 & J_3 \end{bmatrix}, &
 *     {\bf J}_\text{ref} &=
 *         \begin{bmatrix}
 *             J_{\text{ref},0} & J_{\text{ref},1} \\
 *             J_{\text{ref},2} & J_{\text{ref},3}
 *         \end{bmatrix},
 * \end{aligned}
 * \f]
 * then this function computes
 * \f[
 * \begin{bmatrix}
 *     J_0\frac{|J_{\text{ref},0}|}{J_{\text{ref},0}} &
 *     J_1\frac{|J_{\text{ref},1}|}{J_{\text{ref},1}} \\
 *     J_2\frac{|J_{\text{ref},2}|}{J_{\text{ref},2}} &
 *     J_3\frac{|J_{\text{ref},3}|}{J_{\text{ref},3}}.
 * \end{bmatrix}
 * \f]
 *
 * This operation does not affect the magnitudes of the elements of
 * \f${\bf J}\f$, but only their phases.
 */
void remove_reference_phase( gpuDoubleComplex *J, gpuDoubleComplex *Jref )
{
    gpuDoubleComplex PP0norm, PQ0norm, QP0norm, QQ0norm;
    double PPscale = 1.0/gpuCabs( Jref[0] ); // = 1/|PP|
    double PQscale = 1.0/gpuCabs( Jref[1] ); // = 1/|PQ|
    double QPscale = 1.0/gpuCabs( Jref[2] ); // = 1/|QP|
    double QQscale = 1.0/gpuCabs( Jref[3] ); // = 1/|QQ|

    PP0norm = make_gpuDoubleComplex( PPscale*gpuCreal(Jref[0]), PPscale*gpuCimag(Jref[0]) ); // = PP/|PP|
    PQ0norm = make_gpuDoubleComplex( PQscale*gpuCreal(Jref[1]), PQscale*gpuCimag(Jref[1]) ); // = PQ/|PQ|
    QP0norm = make_gpuDoubleComplex( QPscale*gpuCreal(Jref[2]), QPscale*gpuCimag(Jref[2]) ); // = QP/|QP|
    QQ0norm = make_gpuDoubleComplex( QQscale*gpuCreal(Jref[3]), QQscale*gpuCimag(Jref[3]) ); // = QQ/|QQ|

    // Essentially phase rotations
    if (isfinite(PPscale))  J[0] = gpuCdiv( J[0], PP0norm );
    if (isfinite(PQscale))  J[1] = gpuCdiv( J[1], PQ0norm );
    if (isfinite(QPscale))  J[2] = gpuCdiv( J[2], QP0norm );
    if (isfinite(QQscale))  J[3] = gpuCdiv( J[3], QQ0norm );
}

/**
 * Zeroes the off-diagonal terms of the given matrix.
 *
 * @param[in,out] J A complex-valued 2x2 matrix, \f${\bf J}\f$
 *
 * This function sets the off-diagonal terms to zero, thus:
 * \f[
 *     {\bf J} =
 *     \begin{bmatrix} J_0 & J_1 \\ J_2 & J_3 \end{bmatrix}
 *     \rightarrow
 *     \begin{bmatrix} J_0 & 0 \\ 0 & J_3 \end{bmatrix}.
 * \f]
 */
void zero_PQ_and_QP( gpuDoubleComplex *J )
{
    J[1] = make_gpuDoubleComplex( 0.0, 0.0 );
    J[2] = make_gpuDoubleComplex( 0.0, 0.0 );
}

/**
 * Parses the PQ phase correction for the given GPS time from the file
 * `pq_phase_correction.txt`.
 *
 * @param gpstime The GPS second to search for
 * @param cal     The calibration struct to store the read-in information
 */
void parse_calibration_correction_file( uint32_t gpstime, calibration *cal )
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

/**
 * Applies various corrections/adjustments to the calibration solutions.
 *
 * This function applies to the calibration solutions given in `vm&rarr;D`.
 * The corrections to be applied are taken from the calibration struct `vm&rarr;cal`.
 *
 * Three (optional) corrections are applied:
 *   1. Subtract the phases from each antenna by a reference antenna
 *      (see remove_reference_phase()).
 *   2. Remove the PQ and QP (i.e. off-diagonal) terms (see zero_PQ_and_QP()).
 *   3. Apply a phase slope to the QQ terms (implemented in this function).
 */
void vmApplyCalibrationCorrections( vcsbeam_context *vm )
{
    int coarse_chan_idx = vm->coarse_chan_idx;

    // Three locally defined booleans for whether to do the corrections
    bool apply_ref_ant         = true; // Whether this should be true is checked below
    bool apply_zero_PQ_and_QP  = !(vm->cal.keep_cross_terms);
    bool apply_phase_slope     = (vm->cal.phase_offset != 0.0 || vm->cal.phase_slope != 0.0);

    Antenna *Aref = NULL; // The reference antenna struct
    if (strcmp(vm->cal.ref_ant, "NONE") == 0) // i.e. the user explicitly wanted to NOT apply the reference antenna correction
    {
        apply_ref_ant = false;
    }
    else
    {
        Aref = find_antenna_by_name( vm->obs_metadata, vm->cal.ref_ant );

        // If the lookup failed, then quit with an error
        if (Aref == NULL)
        {
            fprintf( stderr, "error: apply_calibration_corrections: "
                    "tile %s does not exist in this observation\n", vm->cal.ref_ant );
            exit(EXIT_FAILURE);
        }
    }

    // Quick check: if all three corrections are turned off, then abort here
    if (!apply_ref_ant && !apply_zero_PQ_and_QP && !apply_phase_slope)
        return;

    // Report on what's being done, if requested
    if (apply_ref_ant | apply_zero_PQ_and_QP | apply_phase_slope)
    {
        sprintf( vm->log_message, "Applying calibration solution corrections to coarse channel %d", coarse_chan_idx );
        logger_timed_message( vm->log, vm->log_message );
    }

    if (apply_ref_ant)
    {
        sprintf( vm->log_message, "    Using tile %s as reference tile for rotating phases", vm->cal.ref_ant );
        logger_timed_message( vm->log, vm->log_message );
    }
    else
    {
        logger_timed_message( vm->log, "    Rotation of phases relative to reference tile is turned OFF" );
    }

    if (apply_zero_PQ_and_QP)
        logger_timed_message( vm->log, "    Setting off-diagonal terms of calibration Jones matrix (PQ and QP) to zero" );
    else
        logger_timed_message( vm->log, "    Retaining off-diagonal terms of calibration Jones matrix (PQ and QP)" );

    if (apply_phase_slope)
    {
        sprintf( vm->log_message, "    Applying phase slope %e*FREQ + %e (rad) to PP", vm->cal.phase_slope, vm->cal.phase_offset );
        logger_timed_message( vm->log, vm->log_message );
    }
    else
        logger_timed_message( vm->log, "    Phase slope correction turned OFF" );

    // Variables for converting slope and offset to a complex phase
    //     z = exp(i*phi) = cos(phi) + i*sin(phi),
    // where
    //     phi = slope*freq + offset

    double          phi; // rad
    gpuDoubleComplex z;   // complex phase

    long int freq_ch; // Hz
    long int frequency  = vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_start_hz; // Hz
    int      chan_width = vm->obs_metadata->coarse_chan_width_hz / vm->nfine_chan;

    uintptr_t d_idx, dref_idx;

    uintptr_t nant    = vm->obs_metadata->num_ants;
    uintptr_t nantpol = vm->obs_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t nchan   = vm->nfine_chan;

    // A temporary copy of the reference antenna matrix, so that we don't clobber it midway
    // through this operation!
    gpuDoubleComplex Dref[nantpol*nantpol];

    // Go through all the antennas and divide the phases by the reference antenna
    uintptr_t ant, ch;
    for (ch = 0; ch < nchan; ch++)
    {
        if (apply_ref_ant)
        {
            // Make a copy of the reference Jones matrix for this channel
            // reference antenna and channel
            dref_idx = D_IDX(Aref->ant,ch,0,0,nchan,nantpol);
            cp2x2( &(vm->D[dref_idx]), Dref );
        }

        if (apply_phase_slope)
        {
            // Convert the slope and offset into a complex phase
            freq_ch = frequency + ch*chan_width;                // The frequency of this fine channel (Hz)
            phi = vm->cal.phase_slope*freq_ch + vm->cal.phase_offset; // (rad)
            z = make_gpuDoubleComplex( cos(phi), sin(phi) );
        }

        for (ant = 0; ant < nant; ant++)
        {
            // A pointer to the (first element of) the Jones matrix for this
            // antenna and channel
            d_idx = D_IDX(ant,ch,0,0,nchan,nantpol);

            // Divide through a reference antenna...
            if (apply_ref_ant)
            {
                remove_reference_phase( &(vm->D[d_idx]), Dref );
            }

            // ...zero the off-diagonal terms...
            if (apply_zero_PQ_and_QP)
                zero_PQ_and_QP( &(vm->D[d_idx]) );

            // ...and apply the phase correction:
            // DZ = [ d00 d01 ] [ z 0 ]
            //      [ d10 d11 ] [ 0 1 ]
            //
            //    = [ d00*z  d01 ]
            //      [ d10*z  d11 ]
            if (apply_phase_slope)
            {
                d_idx = D_IDX(ant,ch,0,0,nchan,nantpol);
                vm->D[d_idx] = gpuCmul( vm->D[d_idx], z );

                d_idx = D_IDX(ant,ch,1,0,nchan,nantpol);
                vm->D[d_idx] = gpuCmul( vm->D[d_idx], z );
            }
        }
    }
}

/**
 * Parse the flagged tile names in the given file.
 *
 * @param filename The name of the file to be parsed
 * @param cal The calibration struct where to store the parsed information
 */
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
    char tilename[MAX_COMMAND_LENGTH]; // More than enough space for tilenames
    cal->nflags = 0;
    while (fscanf( f, "%s", tilename ) != EOF)
        cal->nflags++;

    // Allocate memory for the list of tilenames in the calibration struct
    cal->flagged_tilenames = (char **)malloc( cal->nflags * sizeof(char *) );

    // Rewind to the beginning of the file, and this time read them into
    // the calibration struct
    rewind( f );
    int i;
    for (i = 0; i < cal->nflags; i++)
        fscanf( f, "%ms", &(cal->flagged_tilenames[i]) );

    // Close the file
    fclose( f );
}

/**
 * Tests whether a tile with the given tilename has been flagged.
 *
 * @param tilename The tile name whose flagged status is sought
 * @param cal The calibration struct containing the list of flagged tiles
 *
 * @return True if an only if the tile with the given tilename is flagged in
 *         in `cal`
 */
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

/**
 * Flags tiles by setting the corresponding calibration matrices to zero.
 *
 * If `vm&rarr;cal.flags_file` has been set, this function parses it (by
 * calling vmParseFlaggedTilenamesFile()), and then sets all of the
 * calibration solutions for the named tiles to zero.
 * Afterwards, vmSetNumNotFlaggedRFInputs() is called.
 */
void vmSetCustomTileFlags( vcsbeam_context *vm )
{
    // If no filename given, just count the number of active tiles and return
    if (vm->cal.flags_file == NULL)
    {
        vmSetNumNotFlaggedRFInputs( vm );
        return;
    }

    // Parse the given file for tilenames
    vmParseFlaggedTilenamesFile( vm->cal.flags_file, &vm->cal );

    // Create a "zero" matrix that will be copied
    gpuDoubleComplex Zero[4];
    int i;
    for (i = 0; i < 4; i++)
        Zero[i] = make_gpuDoubleComplex( 0.0, 0.0 );

    // Loop through the tilenames listed in the calibration struct
    char     *tilename; // The tilename in question
    Antenna  *Ant;      // The Antenna struct for a given tilename
    uint32_t  ant;      // The corresponding antenna number
    uintptr_t nchan   = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t nantpol = vm->obs_metadata->num_ant_pols; // = 2 (P, Q)
    uintptr_t ch;       // A particular fine channel
    uintptr_t d_idx;    // Idx into the D array

    for (i = 0; i < vm->cal.nflags; i++)
    {
        // Pointer to the tilename in question
        tilename = vm->cal.flagged_tilenames[i];

        if (!tilename_is_flagged( tilename, &vm->cal ))
            continue;

        // Find the corresponding antenna in the observation
        Ant = find_antenna_by_name( vm->obs_metadata, tilename );
        if (Ant == NULL)
        {
            // No antenna with that name was found, so issue a warning
            // and continue
            sprintf( vm->log_message, "Warning: Tile '%s' (listed in '%s') "
                    "not found\n", tilename, vm->cal.flags_file );
            logger_message( vm->log, vm->log_message );
            continue;
        }
        ant = Ant->ant;

        // Loop through the fine channels
        for (ch = 0; ch < nchan; ch++)
        {
            // Get the index into the D array
            d_idx = D_IDX(ant,ch,0,0,nchan,nantpol);

            // Set it to zero
            cp2x2( Zero, &(vm->D[d_idx]) );
        }
    }

    // Count the number of active tiles
    vmSetNumNotFlaggedRFInputs( vm );
}

/**
 * Sets the member variables of the given calibration struct to their default
 * values.
 *
 * @param cal A pointer to the calibration struct to be initialised
 */
void init_calibration( calibration *cal )
{
    cal->caldir            = NULL;
    cal->ncal_files        = 0;
    cal->cal_files         = NULL;
    cal->ref_ant           = NULL;
    cal->flags_file        = NULL;
    cal->nflags            = 0;
    cal->flagged_tilenames = NULL;
}

/**
 * Frees the memory associated with a calibration struct.
 *
 * @param cal The calibration struct to be freed
 *
 * This function frees memory associated with `cal`, but not `cal` itself.
 */
void free_calibration( calibration *cal )
{
    if (cal->caldir != NULL)
        free( cal->caldir );

    if (cal->cal_files != NULL)
    {
        for (int i = 0; i < cal->ncal_files; i++)
            free( cal->cal_files[i] );
        free( cal->cal_files );
    }

    if (cal->ref_ant != NULL)
        free( cal->ref_ant );

    if (cal->flags_file != NULL)
        free( cal->flags_file );

    if (cal->flagged_tilenames != NULL)
    {
        int i;
        for (i = 0; i < cal->nflags; i++)
        {
            free( cal->flagged_tilenames[i] );
        }
        free( cal->flagged_tilenames );
    }
}
