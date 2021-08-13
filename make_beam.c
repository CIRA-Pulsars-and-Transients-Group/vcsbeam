/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

// TODO: Remove superfluous #includes
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <cuComplex.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <mwalib.h>
#include "ascii_header.h"
#include "mwa_header.h"
#include <glob.h>
#include <fcntl.h>
#include <assert.h>
#include "beam_common.h"
#include "beam_psrfits.h"
#include "beam_vdif.h"
#include "make_beam.h"
#include "vdifio.h"
#include "filter.h"
#include "psrfits.h"
#include "form_beam.h"
#include "calibration.h"
#include <omp.h>

#include <cuda_runtime.h>
#include "ipfb.h"

double now(){
  struct timespec t;
  clock_gettime(CLOCK_REALTIME,&t);
  return (double)t.tv_sec + (double)t.tv_nsec/1000000000L;
}

#define NOW now()

#define ERROR_MESSAGE_LEN  1024

void get_mwalib_metadata(
        struct make_beam_opts *opts,
        MetafitsMetadata **obs_metadata,
        VoltageMetadata  **vcs_metadata,
        VoltageContext   **vcs_context,
        MetafitsMetadata **cal_metadata );

int main(int argc, char **argv)
{
    // A place to hold the beamformer settings
    struct make_beam_opts opts;
    struct calibration cal;           // Variables for calibration settings

    int i; // Generic counter

    /* Set default beamformer settings */

    // Variables for required options
    opts.begin       = 0;    // GPS time -- when to start beamforming
    opts.end         = 0;    // GPS time -- when to stop beamforming
    opts.pointings   = NULL; // list of pointings "dd:mm:ss_hh:mm:ss,dd:mm:ss_hh:mm:ss"
    opts.datadir     = NULL; // The path to where the recombined data live
    opts.metafits    = NULL; // filename of the metafits file for the target observation
    opts.cal_metafits = NULL; // filename of the metafits file for the calibration observation
    opts.rec_channel = -1;   // 0 - 255 receiver 1.28MHz channel

    // Variables for MWA/VCS configuration
    opts.custom_flags  = NULL;   // Use custom list for flagging antennas

    // Output options
    opts.out_incoh     = 0;  // Default = PSRFITS (incoherent) output turned OFF
    opts.out_coh       = 0;  // Default = PSRFITS (coherent)   output turned OFF
    opts.out_vdif      = 0;  // Default = VDIF                 output turned OFF
    opts.out_bf        = 1;  // Default = beamform all (non-flagged) antennas
    opts.out_ant       = 0;  // The antenna number (0-127) to write out if out_bf = 0
    opts.synth_filter  = NULL;
    opts.out_summed    = 0;  // Default = output only Stokes I output turned OFF
    opts.max_sec_per_file = 200; // Number of seconds per fits files

    // Variables for calibration settings
    cal.filename          = NULL;
    cal.bandpass_filename = NULL;
    cal.chan_width        = 40000;
    cal.nchan             = 0;
    cal.cal_type          = NO_CALIBRATION;
    cal.offr_chan_num     = 0;
    cal.ref_ant           = 0;
    cal.cross_terms       = 0;
    cal.phase_offset      = 0.0;
    cal.phase_slope       = 0.0;

    // GPU options
    opts.gpu_mem               = -1.0;

    // Parse command line arguments
    make_beam_parse_cmdline( argc, argv, &opts, &cal );


    double begintime = NOW;

    // Calculate the number of files
    int ntimesteps = opts.end - opts.begin + 1;
    if (ntimesteps <= 0) {
        fprintf(stderr, "Cannot beamform on %d files (between %lu and %lu)\n", ntimesteps, opts.begin, opts.end);
        exit(EXIT_FAILURE);
    }

    // <<<<<
    char error_message[ERROR_MESSAGE_LEN];

    // Create an mwalib metafits context and associated metadata
    fprintf( stderr, "[%f]  Creating metafits and voltage contexts via MWALIB\n", NOW-begintime );

    MetafitsMetadata *obs_metadata = NULL;
    MetafitsMetadata *cal_metadata = NULL;
    VoltageMetadata  *vcs_metadata = NULL;
    VoltageContext   *vcs_context  = NULL;

    if (opts.out_coh || opts.out_vdif)
        get_mwalib_metadata( &opts, &obs_metadata, &vcs_metadata, &vcs_context, &cal_metadata );
    else
        get_mwalib_metadata( &opts, &obs_metadata, &vcs_metadata, &vcs_context, NULL );

    // Create some "shorthand" variables for code brevity
    int nstation             = obs_metadata->num_ants;
    int nchan                = obs_metadata->num_volt_fine_chans_per_coarse;
    int chan_width           = obs_metadata->volt_fine_chan_width_hz;
    int npol                 = obs_metadata->num_ant_pols;   // (X,Y)
    int outpol_coh           = 4;  // (I,Q,U,V)
    if ( opts.out_summed )
        outpol_coh           = 1;  // (I)
    const int outpol_incoh   = 1;  // ("I")
    unsigned int sample_rate = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;

    int ant;
    int offset;
    unsigned int s;
    int ch, pol;


    // =====

    // Read in info from metafits file
    fprintf( stderr, "[%f]  Reading in metafits file information from %s\n", NOW-begintime, opts.metafits);
    struct metafits_info mi;
    get_metafits_info( opts.metafits, &mi, chan_width );

    // >>>>>

    // Start counting time from here (i.e. after parsing the command line)
    fprintf( stderr, "[%f]  Starting %s with GPU acceleration\n", NOW-begintime, argv[0] );

    // Parse input pointings
    int max_npointing = 120; // Could be more
    char RAs[max_npointing][64];
    char DECs[max_npointing][64];
    int npointing = sscanf( opts.pointings,
            "%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,]," ,
                            RAs[0],  DECs[0],  RAs[1],  DECs[1],  RAs[2],  DECs[2],
                            RAs[3],  DECs[3],  RAs[4],  DECs[4],  RAs[5],  DECs[5],
                            RAs[6],  DECs[6],  RAs[7],  DECs[7],  RAs[8],  DECs[8],
                            RAs[9],  DECs[9],  RAs[10], DECs[10], RAs[11], DECs[11],
                            RAs[12], DECs[12], RAs[13], DECs[13], RAs[14], DECs[14],
                            RAs[15], DECs[15], RAs[16], DECs[16], RAs[17], DECs[17],
                            RAs[18], DECs[18], RAs[19], DECs[19], RAs[20], DECs[20],
                            RAs[21], DECs[21], RAs[22], DECs[22], RAs[23], DECs[23],
                            RAs[24], DECs[24], RAs[25], DECs[25], RAs[26], DECs[26],
                            RAs[27], DECs[27], RAs[28], DECs[28], RAs[29], DECs[29],
                            RAs[30], DECs[30], RAs[31], DECs[31], RAs[32], DECs[32],
                            RAs[33], DECs[33], RAs[34], DECs[34], RAs[35], DECs[35],
                            RAs[36], DECs[36], RAs[37], DECs[37], RAs[38], DECs[38],
                            RAs[39], DECs[39], RAs[40], DECs[40], RAs[41], DECs[41],
                            RAs[42], DECs[42], RAs[43], DECs[43], RAs[44], DECs[44],
                            RAs[45], DECs[45], RAs[46], DECs[46], RAs[47], DECs[47],
                            RAs[48], DECs[48], RAs[49], DECs[49], RAs[50], DECs[50],
                            RAs[51], DECs[51], RAs[52], DECs[52], RAs[53], DECs[53],
                            RAs[54], DECs[54], RAs[55], DECs[55], RAs[56], DECs[56],
                            RAs[57], DECs[57], RAs[58], DECs[58], RAs[59], DECs[59],
                            RAs[60], DECs[60], RAs[61], DECs[61], RAs[62], DECs[62],
                            RAs[63], DECs[63], RAs[64], DECs[64], RAs[65], DECs[65],
                            RAs[66], DECs[66], RAs[67], DECs[67], RAs[68], DECs[68],
                            RAs[69], DECs[69], RAs[70], DECs[70], RAs[71], DECs[71],
                            RAs[72], DECs[72], RAs[73], DECs[73], RAs[74], DECs[74],
                            RAs[75], DECs[75], RAs[76], DECs[76], RAs[77], DECs[77],
                            RAs[78], DECs[78], RAs[79], DECs[79], RAs[80], DECs[80],
                            RAs[81], DECs[81], RAs[82], DECs[82], RAs[83], DECs[83],
                            RAs[84], DECs[84], RAs[85], DECs[85], RAs[86], DECs[86],
                            RAs[87], DECs[87], RAs[88], DECs[88], RAs[89], DECs[89],
                            RAs[90], DECs[90], RAs[91], DECs[91], RAs[92], DECs[92],
                            RAs[93], DECs[93], RAs[94], DECs[94], RAs[95], DECs[95],
                            RAs[96], DECs[96], RAs[97], DECs[97], RAs[98], DECs[98],
                            RAs[99], DECs[99], RAs[100], DECs[100], RAs[101], DECs[101],
                            RAs[102], DECs[102], RAs[103], DECs[103], RAs[104], DECs[104],
                            RAs[105], DECs[105], RAs[106], DECs[106], RAs[107], DECs[107],
                            RAs[108], DECs[108], RAs[109], DECs[109], RAs[110], DECs[110],
                            RAs[111], DECs[111], RAs[112], DECs[112], RAs[113], DECs[113],
                            RAs[114], DECs[114], RAs[115], DECs[115], RAs[116], DECs[116],
                            RAs[117], DECs[117], RAs[118], DECs[118], RAs[119], DECs[119] );

    if (npointing%2 == 1)
    {
        fprintf(stderr, "Number of RAs do not equal the number of Decs given. Exiting\n");
        fprintf(stderr, "npointings : %d\n", npointing);
        fprintf(stderr, "RAs[0] : %s\n", RAs[0]);
        fprintf(stderr, "DECs[0] : %s\n", DECs[0]);
        exit(0);
    }
    else
        npointing /= 2; // converting from number of RAs and DECs to number of pointings

    char pointing_array[npointing][2][64];
    int p;
    for ( p = 0; p < npointing; p++)
    {
       strcpy( pointing_array[p][0], RAs[p] );
       strcpy( pointing_array[p][1], DECs[p] );
       fprintf(stderr, "[%f]  Pointing Num: %i  RA: %s  Dec: %s\n", NOW-begintime,
                             p, pointing_array[p][0], pointing_array[p][1]);
    }

    // Allocate memory
    cuDoubleComplex  ****complex_weights_array = create_complex_weights( npointing, nstation, nchan, npol );
    cuDoubleComplex  ****invJi                 = create_invJi( nstation, nchan, npol );
    cuDoubleComplex  ****detected_beam         = create_detected_beam( npointing, 2*sample_rate, nchan, npol );

    // Load the FEE2016 beam and set up the "delays" and "amps" arrays
    fprintf( stderr, "[%f]  Reading in beam model from %s\n", NOW-begintime, HYPERBEAM_HDF5 );
    FEEBeam *beam = NULL;
    beam = new_fee_beam( HYPERBEAM_HDF5 );
    // Also, to prep for the FEE beam function call, create an "amps" array
    // based on the delays array from the metafits
    //int ninputs = obs_metadata->num_rf_inputs;
    int **delays;
    double **amps;
    create_delays_amps_from_metafits( obs_metadata, &delays, &amps );


    // If using bandpass calibration solutions, calculate number of expected bandpass channels
    if (cal.cal_type == RTS_BANDPASS)
        cal.nchan = (nchan * chan_width) / cal.chan_width;

    // If a custom flag file has been provided, use that instead of the metafits flags
    if (opts.custom_flags != NULL)
    {
        // Reset the weights to 1
        for (i = 0; i < nstation*npol; i++)
            mi.weights_array[i] = 1.0;

        // Open custom flag file for reading
        FILE *flagfile = fopen( opts.custom_flags, "r" );
        if (flagfile == NULL)
        {
            fprintf( stderr, "error: couldn't open flag file \"%s\" for "
                             "reading\n", opts.custom_flags );
            exit(EXIT_FAILURE);
        }

        // Read in flags from file
        int nitems;
        int flag;
        while (!feof(flagfile))
        {
            // Read in next item
            nitems = fscanf( flagfile, "%d", &ant );
            if (nitems != 1 && !feof(flagfile))
            {
                fprintf( stderr, "error: couldn't parse flag file \"%s\"\n",
                        opts.custom_flags );
                exit(EXIT_FAILURE);
            }

            // Flag both polarisations of the antenna in question
            flag = ant*2;
            mi.weights_array[flag]   = 0.0;
            mi.weights_array[flag+1] = 0.0;
        }

        // Close file
        fclose( flagfile );
    }

    // Issue warnings if any antennas are being used which are flagged in the metafits file
    for (i = 0; i < nstation*npol; i++)
    {
        if (mi.weights_array[i] != 0.0 &&
            mi.flag_array[i]    != 0.0)
        {
            fprintf( stderr, "warning: antenna %3d, pol %d is included even "
                             "though it is flagged in the metafits file\n",
                             i / npol,
                             i % npol );
        }
    }

    double wgt_sum = 0;
    for (i = 0; i < nstation*npol; i++)
        wgt_sum += mi.weights_array[i];
    double invw = 1.0/wgt_sum;

    // GET CALIBRATION SOLUTION
    // ------------------------
    // 1. Allocate memory
    double amp = 0.0;  // Not actually used yet
    cuDoubleComplex  **M  = (cuDoubleComplex ** ) calloc(nstation, sizeof(cuDoubleComplex * )); // Gain in direction of Calibration
    cuDoubleComplex ***Jf = (cuDoubleComplex ***) calloc(nstation, sizeof(cuDoubleComplex **)); // Fitted bandpass solutions
    cuDoubleComplex Jref[npol*npol];            // Calibration Direction
    cuDoubleComplex invJref[npol*npol];

    int *order = (int *)malloc( nstation*sizeof(int) ); // <-- just for OFFRINGA calibration solutions

    for (ant = 0; ant < nstation; ant++) {
        M[ant]  = (cuDoubleComplex * ) calloc(npol*npol,  sizeof(cuDoubleComplex));
        Jf[ant] = (cuDoubleComplex **) calloc(cal.nchan, sizeof(cuDoubleComplex *));
        for (ch = 0; ch < cal.nchan; ch++) { // Only need as many channels as used in calibration solution
            Jf[ant][ch] = (cuDoubleComplex *) calloc(npol*npol, sizeof(cuDoubleComplex));
        }
    }

    // 2. Read files
    if (cal.cal_type == RTS || cal.cal_type == RTS_BANDPASS)
    {

        read_rts_file(M, Jref, &amp, cal.filename); // Read in the RTS DIJones file
        inv2x2(Jref, invJref);

        if  (cal.cal_type == RTS_BANDPASS) {

            read_bandpass_file(              // Read in the RTS Bandpass file
                    NULL,                    // Output: measured Jones matrices (Jm[ant][ch][pol,pol])
                    Jf,                      // Output: fitted Jones matrices   (Jf[ant][ch][pol,pol])
                    cal.chan_width,         // Input:  channel width of one column in file (in Hz)
                    cal.nchan,              // Input:  (max) number of channels in one file (=128/(chan_width/10000))
                    nstation,                    // Input:  (max) number of antennas in one file (=128)
                    cal.bandpass_filename); // Input:  name of bandpass file

        }

    }
    else if (cal.cal_type == OFFRINGA) {

        // Find the ordering of antennas in Offringa solutions from metafits file
        for (ant = 0; ant < nstation; ant++)
        {
            Antenna A = obs_metadata->antennas[ant];
            order[A.ant*2] = ant;
        }
        read_offringa_gains_file(M, nstation, cal.offr_chan_num, cal.filename, order);
        free(order);

        // Just make Jref (and invJref) the identity matrix since they are already
        // incorporated into Offringa's calibration solutions.
        //Jref[0] = make_cuDoubleComplex( 1.0, 0.0 );
        //Jref[1] = make_cuDoubleComplex( 0.0, 0.0 );
        //Jref[2] = make_cuDoubleComplex( 0.0, 0.0 );
        //Jref[3] = make_cuDoubleComplex( 1.0, 0.0 );
        inv2x2(Jref, invJref);
    }

    // 3. In order to mitigate errors introduced by the calibration scheme, the calibration
    // solution Jones matrix for each antenna may be altered in the following ways:
    //
    //   1) The XX and YY terms are phase-rotated so that those of the supplied
    //      reference antennas are aligned.

    if (cal.ref_ant != -1) // -1 = don't do any phase rotation
        remove_reference_phase( M, cal.ref_ant, nstation );

    //   2) The XY and YX terms are set to zero.

    if (cal.cross_terms == 0)
        zero_XY_and_YX( M, nstation );

    // ------------------------

    // Run get_delays to populate the beam_geom_vals struct
    fprintf( stderr, "[%f]  Setting up output header information\n", NOW-begintime);
    struct beam_geom beam_geom_vals[npointing];

    int coarse_chan_idx = 0; /* Value is fixed for now (i.e. each call of make_beam only
                                ever processes one coarse chan. However, in the future,
                                this should be flexible, with mpi or threads managing
                                different coarse channels. */

    int coarse_chan = vcs_metadata->provided_coarse_chan_indices[coarse_chan_idx];

    get_delays(
            pointing_array,     // an array of pointings [pointing][ra/dec][characters]
            npointing,          // number of pointings
            vcs_metadata,
            obs_metadata,
            coarse_chan_idx,
            &cal,          // struct holding info about calibration
            M,                  // Calibration Jones matrix information
            Jf,                 // Calibration Jones matrix information
            invJref,            // Calibration Jones matrix information
            sample_rate,        // in Hz
            beam,               // Hyperbeam struct
            delays,             // } Analogue beamforming pointing direction information needed for Hyperbeam
            amps,               // }
            0.0,                // seconds offset from the beginning of the observation at which to calculate delays
            beam_geom_vals,     // Populate psrfits header info
            NULL,               // complex weights array (ignore this time)
            NULL                // invJi array           (ignore this time)
    );

    // Create structures for holding header information
    struct psrfits  *pf;
    struct psrfits  *pf_incoh;
    pf = (struct psrfits *)malloc(npointing * sizeof(struct psrfits));
    pf_incoh = (struct psrfits *)malloc(1 * sizeof(struct psrfits));
    vdif_header     vhdr;
    struct vdifinfo *vf;
    vf = (struct vdifinfo *)malloc(npointing * sizeof(struct vdifinfo));


    // Create structures for the PFB filter coefficients
    int ntaps, fil_size = 0;
    double *coeffs = NULL;

    // If no synthesis filter was explicitly chosen, choose the LSQ12 filter
    if (!opts.synth_filter)  opts.synth_filter = strdup("LSQ12");
    if (strcmp( opts.synth_filter, "LSQ12" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchan; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = LSQ12_FILTER_COEFFS; // I'll have to change the way these coefficients are stored
                                                   // in order to avoid this cumbersome loading procedure
        for (i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else if (strcmp( opts.synth_filter, "MIRROR" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchan; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = MIRROR_FILTER_COEFFS;
        for (i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else
    {
        fprintf( stderr, "error: unrecognised synthesis filter: %s\n",
                opts.synth_filter );
        exit(EXIT_FAILURE);
    }
    cuDoubleComplex *twiddles = roots_of_unity( nchan );

    // Adjust by the scaling that was introduced by the forward PFB,
    // along with any other scaling that I, Lord and Master of the inverse
    // PFB, feels is appropriate.
    double approx_filter_scale = 15.0/7.2; // 7.2 = 16384/117964.8
    for (i = 0; i < fil_size; i++)
        coeffs[i] *= approx_filter_scale;

    // Populate the relevant header structs
    populate_psrfits_header( pf,       opts.metafits, obs_metadata->obs_id,
            mi.date_obs, sample_rate, opts.max_sec_per_file, obs_metadata->metafits_coarse_chans[coarse_chan].chan_start_hz, nchan,
            chan_width,outpol_coh, opts.rec_channel, beam_geom_vals,
            mi, npointing, 1 );
    populate_psrfits_header( pf_incoh, opts.metafits, obs_metadata->obs_id,
            mi.date_obs, sample_rate, opts.max_sec_per_file, obs_metadata->metafits_coarse_chans[coarse_chan].chan_start_hz, nchan,
            chan_width, outpol_incoh, opts.rec_channel, beam_geom_vals,
            mi, 1, 0 );

    populate_vdif_header( vf, &vhdr, opts.metafits, obs_metadata->obs_id,
            mi.date_obs, sample_rate, obs_metadata->metafits_coarse_chans[coarse_chan].chan_start_hz, nchan,
            chan_width, opts.rec_channel, beam_geom_vals, npointing );

    // To run asynchronously we require two memory allocations for each data
    // set so multiple parts of the memory can be worked on at once.
    // We control this by changing the pointer to alternate between
    // the two memory allocations

    // Create array for holding the raw data
    long int bytes_per_timestep = vcs_metadata->num_voltage_blocks_per_timestep * vcs_metadata->voltage_block_size_bytes;

    //cudaMallocHost( (void**)&data, bytes_per_timestep * sizeof(uint8_t) );
    uint8_t *data = (uint8_t *)malloc( bytes_per_timestep * sizeof(uint8_t) );
    assert(data);

    /* Allocate host and device memory for the use of the cu_form_beam function */
    // Declaring pointers to the structs so the memory can be alternated
    struct gpu_formbeam_arrays gf;

    struct gpu_ipfb_arrays gi;
    int nchunk;
    malloc_formbeam( &gf, sample_rate, nstation, nchan, npol, &nchunk, opts.gpu_mem,
                     outpol_coh, outpol_incoh, npointing, NOW-begintime );

    // Create a lists of rf_input indexes ordered by antenna number (needed for gpu kernels)
    create_antenna_lists( obs_metadata, gf.polX_idxs, gf.polY_idxs );

    // ... and upload them to the gpu, ready for use!
    cu_upload_pol_idx_lists( &gf );

    // Create output buffer arrays
    float *data_buffer_coh    = NULL;
    float *data_buffer_incoh  = NULL;
    float *data_buffer_vdif   = NULL;

    data_buffer_coh   = create_pinned_data_buffer_psrfits( npointing * nchan *
                                                           outpol_coh * pf[0].hdr.nsblk );
    data_buffer_incoh = create_pinned_data_buffer_psrfits( nchan * outpol_incoh *
                                                           pf_incoh[0].hdr.nsblk );
    data_buffer_vdif  = create_pinned_data_buffer_vdif( vf->sizeof_buffer *
                                                        npointing );

    if (opts.out_vdif)
    {
        malloc_ipfb( &gi, ntaps, sample_rate, nchan, npol, fil_size, npointing );
        cu_load_filter( coeffs, twiddles, &gi, nchan );
    }

    // Set up parrel streams
    cudaStream_t streams[npointing];
    for ( p = 0; p < npointing; p++ )
        cudaStreamCreate(&(streams[p])) ;

    fprintf( stderr, "[%f]  *****BEGINNING BEAMFORMING*****\n", NOW-begintime );

    // Set up timing for each section
    long read_time[ntimesteps], delay_time[ntimesteps], calc_time[ntimesteps], write_time[ntimesteps][npointing];
    int timestep_idx;
    long timestep;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        timestep = vcs_metadata->provided_timestep_indices[timestep_idx];
        coarse_chan = vcs_metadata->provided_coarse_chan_indices[coarse_chan_idx];

        // Read in data from next file
        clock_t start = clock();
        fprintf( stderr, "[%f] [%d/%d] Reading in data for gps second %ld \n", NOW-begintime,
                timestep_idx+1, ntimesteps, vcs_metadata->timesteps[timestep].gps_time_ms / 1000 );

        if (mwalib_voltage_context_read_file(
                    vcs_context,
                    timestep,
                    coarse_chan,
                    data,
                    bytes_per_timestep,
                    error_message,
                    ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
        {
            fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s\n", error_message );
            exit(EXIT_FAILURE);
        }
        read_time[timestep_idx] = clock() - start;

        // Get the next second's worth of phases / jones matrices, if needed
        start = clock();
        fprintf( stderr, "[%f] [%d/%d] Calculating delays\n", NOW-begintime,
                                timestep_idx+1, ntimesteps );
        get_delays(
                pointing_array,     // an array of pointings [pointing][ra/dec][characters]
                npointing,          // number of pointings
                vcs_metadata,
                obs_metadata,
                coarse_chan_idx,
                &cal,              // struct holding info about calibration
                M,                      // Calibration Jones matrix information
                Jf,                     // Calibration Jones matrix information
                invJref,                // Calibration Jones matrix information
                sample_rate,            // Hz
                beam,                   // Hyperbeam struct
                delays,                 // } Analogue beamforming pointing direction information needed for Hyperbeam
                amps,                   // }
                (double)(timestep_idx + opts.begin - obs_metadata->obs_id),        // seconds offset from the beginning of the obseration at which to calculate delays
                NULL,                   // Don't update beam_geom_vals
                complex_weights_array,  // complex weights array (answer will be output here)
                invJi );                // invJi array           (answer will be output here)
        delay_time[timestep_idx] = clock() - start;


        fprintf( stderr, "[%f] [%d/%d] Calculating beam\n", NOW-begintime,
                                timestep_idx+1, ntimesteps);
        start = clock();

        if (!opts.out_bf) // don't beamform, but only procoess one ant/pol combination
        {
            // Populate the detected_beam, data_buffer_coh, and data_buffer_incoh arrays
            // detected_beam = [2*sample_rate][nchan][npol] = [20000][128][2] (cuDoubleComplex)
            if (timestep_idx % 2 == 0)
                offset = 0;
            else
                offset = sample_rate;

            for (p   = 0; p   < npointing;   p++  )
            for (s   = 0; s   < sample_rate; s++  )
            for (ch  = 0; ch  < nchan;       ch++ )
            for (pol = 0; pol < npol;        pol++)
            {
                if (pol == 0)
                    i = gf.polX_idxs[opts.out_ant];
                else
                    i = gf.polY_idxs[opts.out_ant];
                detected_beam[p][s+offset][ch][pol] = UCMPLX4_TO_CMPLX_FLT(data[D_IDX(s,ch,i,nchan)]);
                detected_beam[p][s+offset][ch][pol] = cuCmul(detected_beam[p][s+offset][ch][pol], make_cuDoubleComplex(128.0, 0.0));
            }
        }
        else // beamform (the default mode)
        {
            cu_form_beam( data, sample_rate, complex_weights_array, invJi, timestep_idx,
                    npointing, nstation, nchan, npol, outpol_coh, invw, &gf,
                    detected_beam, data_buffer_coh, data_buffer_incoh,
                    streams, opts.out_incoh, nchunk );
        }

        // Invert the PFB, if requested
        if (opts.out_vdif)
        {
            fprintf( stderr, "[%f] [%d/%d] Inverting the PFB (full)\n",
                            NOW-begintime, timestep_idx+1, ntimesteps);
            cu_invert_pfb_ord( detected_beam, timestep_idx, npointing,
                               sample_rate, nchan, npol, vf->sizeof_buffer,
                               &gi, data_buffer_vdif );
        }
        calc_time[timestep_idx] = clock() - start;


        // Write out for each pointing
        for ( p = 0; p < npointing; p++)
        {
            start = clock();
            fprintf( stderr, "[%f] [%d/%d] [%d/%d] Writing data to file(s)\n",
                    NOW-begintime, timestep_idx+1, ntimesteps, p+1, npointing );

            if (opts.out_coh)
                psrfits_write_second( &pf[p], data_buffer_coh, nchan,
                                      outpol_coh, p );
            if (opts.out_incoh && p == 0)
                psrfits_write_second( &pf_incoh[p], data_buffer_incoh,
                                      nchan, outpol_incoh, p );
            if (opts.out_vdif)
                vdif_write_second( &vf[p], &vhdr,
                                   data_buffer_vdif + p * vf->sizeof_buffer );
            write_time[timestep_idx][p] = clock() - start;
        }
    }

    // Calculate total processing times
    float read_sum = 0, delay_sum = 0, calc_sum = 0, write_sum = 0;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        read_sum  += (float) read_time[timestep_idx];
        delay_sum += (float) delay_time[timestep_idx];
        calc_sum  += (float) calc_time[timestep_idx];
        for ( p = 0; p < npointing; p++)
            write_sum += (float) write_time[timestep_idx][p];
    }
    float read_mean, delay_mean, calc_mean, write_mean;
    read_mean  = read_sum  / ntimesteps;
    delay_mean = delay_sum / ntimesteps / npointing;
    calc_mean  = calc_sum  / ntimesteps / npointing;
    write_mean = write_sum / ntimesteps / npointing;

    // Calculate the standard deviations
    float read_std = 0, delay_std = 0, calc_std = 0, write_std = 0;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        read_std  += pow((float)read_time[timestep_idx]  - read_mean,  2);
        delay_std += pow((float)delay_time[timestep_idx] - delay_mean / npointing, 2);
        calc_std  += pow((float)calc_time[timestep_idx]  - calc_mean / npointing,  2);
        for ( p = 0; p < npointing; p++)
            write_std += pow((float)write_time[timestep_idx][p] - write_mean / npointing, 2);
    }
    read_std  = sqrt( read_std  / ntimesteps );
    delay_std = sqrt( delay_std / ntimesteps / npointing );
    calc_std  = sqrt( calc_std  / ntimesteps / npointing );
    write_std = sqrt( write_std / ntimesteps / npointing );


    fprintf( stderr, "[%f]  **FINISHED BEAMFORMING**\n", NOW-begintime);
    fprintf( stderr, "[%f]  Total read  processing time: %9.3f s\n",
                     NOW-begintime, read_sum / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  read  processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, read_mean / CLOCKS_PER_SEC, read_std / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Total delay processing time: %9.3f s\n",
                     NOW-begintime, delay_sum / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  delay processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, delay_mean / CLOCKS_PER_SEC, delay_std / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Total calc  processing time: %9.3f s\n",
                     NOW-begintime, calc_sum / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  calc  processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, calc_mean / CLOCKS_PER_SEC, calc_std / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Total write processing time: %9.3f s\n",
                     NOW-begintime, write_sum  * npointing / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  write processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, write_mean / CLOCKS_PER_SEC, write_std / CLOCKS_PER_SEC);


    fprintf( stderr, "[%f]  Starting clean-up\n", NOW-begintime);

    // Free up memory
    destroy_complex_weights( complex_weights_array, npointing, nstation, nchan );
    destroy_invJi( invJi, nstation, nchan, npol );
    destroy_detected_beam( detected_beam, npointing, 2*sample_rate, nchan );

    free( twiddles );
    free( coeffs );

    destroy_metafits_info( &mi );
    //free( data_buffer_coh    );
    //free( data_buffer_incoh  );
    //free( data_buffer_vdif   );
    cudaFreeHost( data_buffer_coh   );
    cudaFreeHost( data_buffer_incoh );
    cudaFreeHost( data_buffer_vdif  );
    cudaFreeHost( data );

    free( opts.pointings    );
    free( opts.datadir      );
    free( opts.metafits     );
    free( cal.filename );
    free( opts.custom_flags );
    free( opts.synth_filter );

    if (opts.out_incoh)
    {
        free( pf_incoh[0].sub.data        );
        free( pf_incoh[0].sub.dat_freqs   );
        free( pf_incoh[0].sub.dat_weights );
        free( pf_incoh[0].sub.dat_offsets );
        free( pf_incoh[0].sub.dat_scales  );
    }
    for (p = 0; p < npointing; p++)
    {
        if (opts.out_coh)
        {
            free( pf[p].sub.data        );
            free( pf[p].sub.dat_freqs   );
            free( pf[p].sub.dat_weights );
            free( pf[p].sub.dat_offsets );
            free( pf[p].sub.dat_scales  );
        }
        if (opts.out_vdif)
        {
            free( vf[p].b_scales  );
            free( vf[p].b_offsets );
        }
    }
    free_formbeam( &gf );
    if (opts.out_vdif)
    {
        free_ipfb( &gi );
    }

    // Clean up Hyperbeam
    free_fee_beam( beam );
    free_delays_amps( obs_metadata, delays, amps );

    mwalib_metafits_metadata_free( obs_metadata );
    mwalib_voltage_metadata_free( vcs_metadata );
    mwalib_voltage_context_free( vcs_context );
    if (cal_metadata != NULL)
        mwalib_metafits_metadata_free( cal_metadata );

    // Clean up memory used for calibration solutions
    for (ant = 0; ant < nstation; ant++) {
        for (ch = 0; ch < cal.nchan; ch++)
            free( Jf[ant][ch] );
        free( Jf[ant] );
        free( M[ant] );
    }

    free( Jf );
    free( M );
    free( order );

    return EXIT_SUCCESS;
}


void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: make_beam [OPTIONS]\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "REQUIRED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-b, --begin=GPSTIME       ");
    fprintf(stderr, "Begin time of observation, in GPS seconds\n");
    fprintf(stderr, "\t-e, --end=GPSTIME         ");
    fprintf(stderr, "End time of observation, in GPS seconds\n");
    fprintf(stderr, "\t-P, --pointings=hh:mm:ss.s_dd:mm:ss.s,hh:mm:ss.s_dd:mm:ss.s...\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "Right ascension and declinations of multiple pointings\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-d, --data-location=PATH  ");
    fprintf(stderr, "PATH is the directory containing the recombined data\n");
    fprintf(stderr, "\t-m, --metafits=FILE  ");
    fprintf(stderr, "FILE is the metafits file for the target observation\n");
    fprintf(stderr, "\t-f, --coarse-chan=N       ");
    fprintf(stderr, "Absolute coarse channel number (0-255)\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-i, --incoh                ");
    fprintf(stderr, "Turn on incoherent PSRFITS beam output.                             ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-p, --psrfits              ");
    fprintf(stderr, "Turn on coherent PSRFITS output (will be turned on if none of\n");
    fprintf(stderr, "\t                           ");
    fprintf(stderr, "-i, -p, -u, -v are chosen).                                         ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-v, --vdif                 ");
    fprintf(stderr, "Turn on VDIF output with upsampling                                 ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-s, --summed               ");
    fprintf(stderr, "Turn on summed polarisations of the coherent output (only Stokes I) ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-t, --max_t                ");
    fprintf(stderr, "Maximum number of seconds per output fits file. ");
    fprintf(stderr, "[default: 200]\n");
    fprintf(stderr, "\t-A, --antpol=ant           ");
    fprintf(stderr, "Do not beamform. Instead, only operate on the specified ant\n");
    fprintf(stderr, "\t                           ");
    fprintf(stderr, "stream (0-127)\n" );
    fprintf(stderr, "\t-S, --synth_filter=filter  ");
    fprintf(stderr, "Apply the named filter during high-time resolution synthesis.    ");
    fprintf(stderr, "[default: LSQ12]\n");
    fprintf(stderr, "\t                           ");
    fprintf(stderr, "filter can be MIRROR or LSQ12.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "MWA/VCS CONFIGURATION OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-F, --custom-flags=file   ");
    fprintf(stderr, "Flag the antennas listed in file instead of those flagged in the ");
    fprintf(stderr, "[default: none]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "metafits file given by the -m option.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "CALIBRATION OPTIONS (RTS)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-c, --cal-metafits=FILE  ");
    fprintf(stderr, "FILE is the metafits file pertaining to the calibration solution\n");
    fprintf(stderr, "\t-J, --dijones-file=PATH   ");
    fprintf(stderr, "The direction-independent Jones matrix file that is output from\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the RTS. Using this option instructs the beamformer to use the\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "RTS-generated calibration solution. Either -J or -O must be\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "supplied. If both are supplied the one that comes last will\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "override the former.\n");
    fprintf(stderr, "\t-B, --bandpass-file=PATH  ");
    fprintf(stderr, "The bandpass file that is output from the RTS. If this option\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "is given, the RTS calibration solution will be applied to each\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "fine channel. If -J is supplied but -B is not, then the coarse\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "channel solution will be applied to ALL fine channels\n");
    fprintf(stderr, "\t-R, --ref-ant=ANT         ");
    fprintf(stderr, "Rotate the phases of the XX and YY elements of the calibration\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "Jones matrices so that the phases of tile ANT align. If ANT is\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "outside the range 0-127, no phase rotation is done               ");
    fprintf(stderr, "[default: 0]\n");
    fprintf(stderr, "\t-X, --cross-terms         ");
    fprintf(stderr, "Retain the XY and YX terms of the calibration solution           ");
    fprintf(stderr, "[default: off]\n");
    fprintf(stderr, "\t-U, --UV-phase=M,C        ");
    fprintf(stderr, "Rotate the Y pol by M*f+C, where M is in rad/Hz and C is in rad  ");
    fprintf(stderr, "[default: 0.0,0.0]\n");
    fprintf(stderr, "\t-W, --rts-chan-width      ");
    fprintf(stderr, "RTS calibration channel bandwidth (Hz)                           ");
    fprintf(stderr, "[default: 40000]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CALIBRATION OPTIONS (OFFRINGA)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-O, --offringa-file=PATH  ");
    fprintf(stderr, "The calibration solution file that is output from the tools\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "made by Andre Offringa. Using this option instructs the beam-\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "former to use the Offringa-style calibration solution. Either\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "-J or -O must be supplied. If both are supplied the one that\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "comes last will override the former.\n");
    fprintf(stderr, "\t-C, --offringa-chan=N     ");
    fprintf(stderr, "The zero-offset position of the coarse channel solution in the   ");
    fprintf(stderr, "[default: 0]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "calibration file given by the -O option.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OTHER OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-h, --help                ");
    fprintf(stderr, "Print this help and exit\n");
    fprintf(stderr, "\t-g, --gpu-mem=N     ");
    fprintf(stderr, "The maximum amount of GPU memory you want make_beam to use in GB ");
    fprintf(stderr, "[default: -1]\n");
    fprintf(stderr, "\t-V, --version             ");
    fprintf(stderr, "Print version number and exit\n");
    fprintf(stderr, "\n");
}



void make_beam_parse_cmdline(
        int argc, char **argv, struct make_beam_opts *opts, struct calibration *cal )
{
    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"end",             required_argument, 0, 'e'},
                {"incoh",           no_argument,       0, 'i'},
                {"psrfits",         no_argument,       0, 'p'},
                {"vdif",            no_argument,       0, 'v'},
                {"summed",          no_argument,       0, 's'},
                {"max_t",           required_argument, 0, 't'},
                {"synth_filter",    required_argument, 0, 'S'},
                {"antpol",          required_argument, 0, 'A'},
                {"pointings",       required_argument, 0, 'P'},
                {"data-location",   required_argument, 0, 'd'},
                {"metafits",        required_argument, 0, 'm'},
                {"cal-metafits",    required_argument, 0, 'c'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"custom-flags",    required_argument, 0, 'F'},
                {"dijones-file",    required_argument, 0, 'J'},
                {"bandpass-file",   required_argument, 0, 'B'},
                {"ref-ant",         required_argument, 0, 'R'},
                {"cross-terms",     no_argument,       0, 'X'},
                {"UV-phase",        required_argument, 0, 'U'},
                {"rts-chan-width",  required_argument, 0, 'W'},
                {"offringa-file",   required_argument, 0, 'O'},
                {"offringa-chan",   required_argument, 0, 'C'},
                {"gpu-mem",         required_argument, 0, 'g'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "A:b:B:c:C:d:e:f:F:g:hiJ:m:O:pP:R:sS:t:U:vVW:X",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'A':
                    opts->out_bf = 0; // Turn off normal beamforming
                    opts->out_ant = atoi(optarg); // 0-127
                    break;
                case 'b':
                    opts->begin = atol(optarg);
                    break;
                case 'B':
                    cal->bandpass_filename = strdup(optarg);
                    cal->cal_type = RTS_BANDPASS;
                    break;
                case 'c':
                    opts->cal_metafits = strdup(optarg);
                    break;
                case 'C':
                    cal->offr_chan_num = atoi(optarg);
                    break;
                case 'd':
                    opts->datadir = strdup(optarg);
                    break;
                case 'e':
                    opts->end = atol(optarg);
                    break;
                case 'f':
                    opts->rec_channel = atoi(optarg);
                    break;
                case 'F':
                    opts->custom_flags = strdup(optarg);
                    break;
                case 'g':
                    opts->gpu_mem = atof(optarg);
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'i':
                    opts->out_incoh = 1;
                    break;
                case 'J':
                    cal->filename = strdup(optarg);
                    if (cal->cal_type != RTS_BANDPASS)
                        cal->cal_type = RTS;
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'O':
                    cal->filename = strdup(optarg);
                    cal->cal_type = OFFRINGA;
                    break;
                case 'p':
                    opts->out_coh = 1;
                    break;
                case 'P':
                    opts->pointings = strdup(optarg);
                    break;
                case 'R':
                    cal->ref_ant = atoi(optarg);
                    break;
                case 'S':
                    opts->synth_filter = strdup(optarg);
                    break;
                case 's':
                    opts->out_summed = 1;
                    break;
                case 't':
                    opts->max_sec_per_file = atoi(optarg);
                    break;
                case 'U':
                    if (sscanf( optarg, "%lf,%lf", &(cal->phase_slope),
                                &(cal->phase_offset) ) != 2)
                    {
                        fprintf( stderr, "error: badly formed argument to -Y "
                                "('%s')\n", optarg );
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'v':
                    opts->out_vdif = 1;
                    break;
                case 'V':
                    fprintf( stderr, "MWA Beamformer %s\n", VERSION_BEAMFORMER);
                    exit(0);
                    break;
                case 'W':
                    cal->chan_width = atoi(optarg);
                    break;
                case 'X':
                    cal->cross_terms = 1;
                    break;
                default:
                    fprintf(stderr, "error: make_beam_parse_cmdline: "
                                    "unrecognised option '%s'\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
            }
        }
    }
    else {
        usage();
        exit(EXIT_FAILURE);
    }


    // Check that all the required options were supplied
    assert( opts->begin        != 0    );
    assert( opts->end          != 0    );
    assert( opts->pointings    != NULL );
    assert( opts->datadir      != NULL );
    assert( opts->metafits     != NULL );
    assert( opts->rec_channel  != -1   );
    assert( cal->cal_type != NO_CALIBRATION );
    if (opts->out_coh || opts->out_vdif)
        assert( opts->cal_metafits );

    // If neither -i, -p, nor -v were chosen, set -p by default
    if ( !opts->out_incoh && !opts->out_coh && !opts->out_vdif )
    {
        opts->out_coh = 1;
    }

    // If the reference antenna is outside the range of antennas, issue
    // a warning and turn off phase rotation.
    /*if (cal->ref_ant < 0 || cal->ref_ant >= opts->nstation)
    {
        fprintf( stderr, "warning: tile %d outside of range 0-%d. "
                "Calibration phase rotation turned off.\n",
                cal->ref_ant, opts->nstation-1 );
        cal->ref_ant = -1; // This defines the meaning of -1
                                // = don't do phase rotation
    }*/
}



char **create_filenames( const struct MetafitsContext *metafits_context, const struct MetafitsMetadata *metafits_metadata, struct make_beam_opts *opts )
// Create an array of filenames; free with destroy_filenames()
{
    // Buffer for mwalib error messages
    char error_message[ERROR_MESSAGE_LEN];

    // Calculate the number of files
    unsigned int ntimesteps = opts->end - opts->begin + 1;
    if (ntimesteps <= 0) {
        fprintf( stderr, "Cannot beamform on %d files (between %ld and %ld)\n",
                 ntimesteps, opts->begin, opts->end);
        exit(EXIT_FAILURE);
    }
    // Allocate memory for the file name list
    char filename[MAX_COMMAND_LENGTH]; // Just the mwalib-generated filename (without the path)
    char **filenames = (char **)malloc( ntimesteps*sizeof(char *) ); // The full array of filenames, including the paths

    // Get the coarse channel index
    unsigned int coarse_chan_idx;
    for (coarse_chan_idx = 0; coarse_chan_idx < metafits_metadata->num_metafits_coarse_chans; coarse_chan_idx++)
        if (metafits_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number == (uintptr_t)opts->rec_channel)
            break;

    // Allocate memory and write filenames
    unsigned int timestep, timestep_idx, filename_idx;
    for (timestep = opts->begin; timestep <= opts->end; timestep++)
    {
        timestep_idx = timestep - (metafits_metadata->metafits_timesteps[0].gps_time_ms / 1000); // <-- This assumes timesteps are always contiguous
        filename_idx = timestep - opts->begin;

        filenames[filename_idx] = (char *)malloc( MAX_COMMAND_LENGTH*sizeof(char) );
        if (mwalib_metafits_get_expected_volt_filename(
                    metafits_context,
                    timestep_idx,
                    coarse_chan_idx,
                    filename,
                    MAX_COMMAND_LENGTH,
                    error_message,
                    ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
        {
            fprintf( stderr, "error (mwalib): cannot create filenames: %s\n", error_message );
            exit(EXIT_FAILURE);
        }
        sprintf( filenames[filename_idx], "%s/%s",
                 opts->datadir, filename );
    }

    return filenames;
}

void destroy_filenames( char **filenames, int nfiles )
{
    int i;
    for (i = 0; i < nfiles; i++)
        free( filenames[i] );
    free( filenames );
}


cuDoubleComplex ****create_complex_weights( int npointing, int nstation, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, ant, ch; // Loop variables
    cuDoubleComplex ****array;

    array = (cuDoubleComplex ****)malloc( npointing * sizeof(cuDoubleComplex ***) );

    for (p = 0; p < npointing; p++)
    {
        array[p] = (cuDoubleComplex ***)malloc( nstation * sizeof(cuDoubleComplex **) );

        for (ant = 0; ant < nstation; ant++)
        {
            array[p][ant] = (cuDoubleComplex **)malloc( nchan * sizeof(cuDoubleComplex *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][ant][ch] = (cuDoubleComplex *)malloc( npol * sizeof(cuDoubleComplex) );
        }
    }
    return array;
}


void destroy_complex_weights( cuDoubleComplex ****array, int npointing, int nstation, int nchan )
{
    int p, ant, ch;
    for (p = 0; p < npointing; p++)
    {
        for (ant = 0; ant < nstation; ant++)
        {
            for (ch = 0; ch < nchan; ch++)
                free( array[p][ant][ch] );

            free( array[p][ant] );
        }
        free( array[p] );
    }
    free( array );
}

cuDoubleComplex ****create_invJi( int nstation, int nchan, int npol )
// Allocate memory for (inverse) Jones matrices
{
    int ant, pol, ch; // Loop variables
    cuDoubleComplex ****invJi;
    invJi = (cuDoubleComplex ****)malloc( nstation * sizeof(cuDoubleComplex ***) );

    for (ant = 0; ant < nstation; ant++)
    {
        invJi[ant] =(cuDoubleComplex ***)malloc( nchan * sizeof(cuDoubleComplex **) );

        for (ch = 0; ch < nchan; ch++)
        {
            invJi[ant][ch] = (cuDoubleComplex **)malloc( npol * sizeof(cuDoubleComplex *) );

            for (pol = 0; pol < npol; pol++)
                invJi[ant][ch][pol] = (cuDoubleComplex *)malloc( npol * sizeof(cuDoubleComplex) );
        }
    }
    return invJi;
}


void destroy_invJi( cuDoubleComplex ****array, int nstation, int nchan, int npol )
{
    int ant, ch, pol;
    for (ant = 0; ant < nstation; ant++)
    {
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
                free( array[ant][ch][pol] );

            free( array[ant][ch] );
        }
        free( array[ant] );
    }
    free( array );
}


cuDoubleComplex ****create_detected_beam( int npointing, int nsamples, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, s, ch; // Loop variables
    cuDoubleComplex ****array;

    array = (cuDoubleComplex ****)malloc( npointing * sizeof(cuDoubleComplex ***) );
    for (p = 0; p < npointing; p++)
    {
        array[p] = (cuDoubleComplex ***)malloc( nsamples * sizeof(cuDoubleComplex **) );

        for (s = 0; s < nsamples; s++)
        {
            array[p][s] = (cuDoubleComplex **)malloc( nchan * sizeof(cuDoubleComplex *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][s][ch] = (cuDoubleComplex *)malloc( npol * sizeof(cuDoubleComplex) );
        }
    }
    return array;
}

void destroy_detected_beam( cuDoubleComplex ****array, int npointing, int nsamples, int nchan )
{
    int p, s, ch;
    for (p = 0; p < npointing; p++)
    {
        for (s = 0; s < nsamples; s++)
        {
            for (ch = 0; ch < nchan; ch++)
                free( array[p][s][ch] );

            free( array[p][s] );
        }

        free( array[p] );
    }

    free( array );
}

float *create_data_buffer_psrfits( size_t size )
{
    float *ptr = (float *)malloc( size * sizeof(float) );
    return ptr;
}


float *create_data_buffer_vdif( size_t size )
{
    float *ptr  = (float *)malloc( size * sizeof(float) );
    return ptr;
}


void get_mwalib_metadata(
        struct make_beam_opts *opts,
        MetafitsMetadata **obs_metadata,
        VoltageMetadata  **vcs_metadata,
        VoltageContext   **vcs_context,
        MetafitsMetadata **cal_metadata )
/* Create the metadata structs using mwalib's API.
 *
 * obs_metadata applies to the metafits metadata for the target observation, and cannot be null
 * vcs_metadata applies to the voltage metadata for the target observation, and cannot be null
 * cal_metadata applies to the metafits metadata for the calibration observation, and may be NULL
 * vcs_context applies to the voltage context for the target observation, and cannot be null
 *     (if only an incoherent beam is requested)
 */
{
    char error_message[ERROR_MESSAGE_LEN];

    // Each metadata is constructed from mwalib "contexts"
    // These are only temporarily needed to construct the other metadata/context structs,
    // and will be freed at the end of this function
    MetafitsContext *obs_context = NULL;
    MetafitsContext *cal_context = NULL;

    // First, get the metafits context for the given observation metafits file in order to create
    // a list of filenames for creating the voltage context. This metafits context will be remade
    // from the voltage context later, in order to get the antenna ordering correct (which depends
    // on the voltage type)

    // Create OBS_CONTEXT
    if (mwalib_metafits_context_new2( opts->metafits, &obs_context, error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( obs_context, NULL, NULL, obs_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create list of filenames
    char **filenames = create_filenames( obs_context, *obs_metadata, opts );

    // Now that we have the filenames, we don't need this version of the obs context and metadata.
    mwalib_metafits_metadata_free( *obs_metadata );
    mwalib_metafits_context_free( obs_context );

    // Create an mwalib voltage context, voltage metadata, and new obs metadata (now with correct antenna ordering)
    // (MWALIB is expecting a const array, so we will give it one!)
    unsigned int ntimesteps = opts->end - opts->begin + 1;
    const char **voltage_files = (const char **)malloc( sizeof(char *) * ntimesteps );
    unsigned int i;
    for (i = 0; i < ntimesteps; i++)
        voltage_files[i] = filenames[i];

    // Create VCS_CONTEXT
    if (mwalib_voltage_context_new( opts->metafits, voltage_files, ntimesteps, vcs_context, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( NULL, NULL, *vcs_context, obs_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata from voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create VCS_METADATA
    if (mwalib_voltage_metadata_get( *vcs_context, vcs_metadata, error_message, ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot get metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Free memory
    destroy_filenames( filenames, ntimesteps );
    free( voltage_files );

    // Can stop at this point if no calibration metadata is requested
    if (cal_metadata == NULL)
        return;

    // Create CAL_CONTEXT
    if (mwalib_metafits_context_new2( opts->cal_metafits, &cal_context, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create cal metafits context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create CAL_METADATA
    if (mwalib_metafits_metadata_get( cal_context, NULL, NULL, cal_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create cal metafits metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }
}
