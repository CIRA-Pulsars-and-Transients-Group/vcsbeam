/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

// Standard library
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>

// Non-standard dependencies
#include <mwalib.h>
#include <mpi.h>

// Local includes
#include "vcsbeam.h"

#define MAX_COMMAND_LENGTH 1024

struct make_tied_array_beam_opts {
    char              *begin_str;        // Absolute or relative GPS time -- when to start beamforming
    unsigned long int  nseconds;         // How many seconds to process
    char              *pointings_file;   // Name of file containing pointings (e.g. "hh:mm:ss dd:mm:ss")
    char              *datadir;          // The path to where the recombined data live
    char              *metafits;         // filename of the metafits file
    char              *coarse_chan_str;  // Absolute or relative coarse channel number

    // Variables for MWA/VCS configuration
    //char              *custom_flags;    // Use custom list for flagging antennas

    // Output options
    bool               out_fine;         // Output fine channelised data (PSRFITS)
    bool               out_coarse;       // Output coarse channelised data (VDIF)
    bool               emulate_legacy;   // Use the offline forward pfb if input data is MWAX

    // Other options
    char              *synth_filter;     // Which synthesis filter to use
    int                max_sec_per_file; // Number of seconds per fits files
    float              gpu_mem;          // Default = -1.0. If -1.0 use all GPU mem
};

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/

void usage();
void make_tied_array_beam_parse_cmdline( int argc, char **argv, struct make_tied_array_beam_opts *opts, struct calibration *cal );
void parse_pointing_file( const char *filename, double **ras_hours, double **decs_degs, unsigned int *npointings );

void write_step( mpi_psrfits *mpfs, int npointing, bool out_fine, bool out_coarse,
        struct vdifinfo *vf, vdif_header *vhdr, float *data_buffer_vdif, logger *log );

/********
 * MAIN *
 ********/

int main(int argc, char **argv)
{
    // Initialise MPI
    MPI_Init( NULL, NULL );
    int world_size, mpi_proc_id;
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_proc_id );

    const int writer = 0; // Designate process 0 to write out the files
    const int ncoarse_chans = world_size;

    // Parse command line arguments
    struct make_tied_array_beam_opts opts;
    struct calibration cal;           // Variables for calibration settings
    make_tied_array_beam_parse_cmdline( argc, argv, &opts, &cal );

    int i; // Generic counter

    // Start a logger for output messages and time-keeping
    logger *log = create_logger( stdout, mpi_proc_id );
    logger_add_stopwatch( log, "read", "Reading in data" );
    logger_add_stopwatch( log, "delay", "Calculating geometric and cable delays" );
    logger_add_stopwatch( log, "calc", "Calculating tied-array beam" );
    logger_add_stopwatch( log, "ipfb", "Inverting the PFB" );
    logger_add_stopwatch( log, "splice", "Splicing coarse channels together" );
    logger_add_stopwatch( log, "write", "Writing out data to file" );
    char log_message[MAX_COMMAND_LENGTH];

    // Create an mwalib metafits context and associated metadata
    logger_timed_message( log, "Creating metafits and voltage contexts via MWALIB" );

    const int chans_per_proc = 1;
    vcsbeam_metadata *vm = init_vcsbeam_metadata(
        opts.metafits, cal.metafits,
        opts.coarse_chan_str, chans_per_proc, mpi_proc_id,
        opts.begin_str, opts.nseconds, 0,
        opts.datadir );

    if (opts.out_fine)    set_vcsbeam_fine_output( vm, true );
    if (opts.out_coarse)  set_vcsbeam_coarse_output( vm, true );

    // Create some "shorthand" variables for code brevity
    uintptr_t nants          = vm->obs_metadata->num_ants;
    uintptr_t nchans         = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t npols          = vm->obs_metadata->num_ant_pols;   // (X,Y)
    unsigned int nsamples    = vm->sample_rate;
    uintptr_t data_size      = vm->bytes_per_second;


    uintptr_t timestep_idx;
    uint64_t  gps_second;

    primary_beam pb;
    geometric_delays gdelays;

    // Start counting time from here (i.e. after parsing the command line)
    sprintf( log_message, "Reading pointings file %s", opts.pointings_file );
    logger_timed_message( log, log_message );

    // Parse input pointings
    double *ras_hours, *decs_degs;
    unsigned int npointing;
    unsigned int p;
    parse_pointing_file( opts.pointings_file, &ras_hours, &decs_degs, &npointing );

    // Allocate memory for various data products
    cuDoubleComplex  ****detected_beam = create_detected_beam( npointing, 2*nsamples, nchans, npols );

    double invw = 1.0/get_num_not_flagged_rf_inputs( vm );

    // Get pointing geometry information
    struct beam_geom beam_geom_vals[npointing];

    double mjd, sec_offset;
    mjd = vm->obs_metadata->sched_start_mjd;
    for (p = 0; p < npointing; p++)
        calc_beam_geom( ras_hours[p], decs_degs[p], mjd, &beam_geom_vals[p] );

    // Create a structure for the PFB filter coefficients

    // If no synthesis filter was explicitly chosen, choose the LSQ12 filter
    if (!opts.synth_filter)
    {
        opts.synth_filter = (char *)malloc( 6 );
        strcpy( opts.synth_filter, "LSQ12" );
    }

    pfb_filter *filter = load_filter_coefficients( opts.synth_filter, SYNTHESIS_FILTER, nchans );

    // Adjust by the scaling that was introduced by the forward PFB,
    // along with any other scaling that I, Lord and Master of the inverse
    // PFB, feel is appropriate.
    double approx_filter_scale = 15.0/7.2; // 7.2 = 16384/117964.8
    for (i = 0; i < filter->ncoeffs; i++)
        filter->coeffs[i] *= approx_filter_scale;

    /*********************
     * Memory Allocation *
     *********************/

    uint8_t *data;
    cudaMallocHost( (void **)&data, data_size );
    cudaCheckErrors( "cudaMallocHost(data) failed" );

    /* Allocate host and device memory for the use of the cu_form_beam function */
    // Declaring pointers to the structs so the memory can be alternated
    struct gpu_formbeam_arrays gf;

    struct gpu_ipfb_arrays gi;
    int nchunk;
    malloc_formbeam( &gf, vm, &nchunk, opts.gpu_mem,
                     NSTOKES, npointing, log );

    // Create a lists of rf_input indexes ordered by antenna number (needed for gpu kernels)
    create_antenna_lists( vm->obs_metadata, gf.polX_idxs, gf.polY_idxs );

    // ... and upload them to the gpu, ready for use!
    cu_upload_pol_idx_lists( &gf );

    // Create output buffer arrays
    float *data_buffer_coh    = NULL;
    float *data_buffer_vdif   = NULL;

    data_buffer_coh   = create_pinned_data_buffer_psrfits( npointing * nchans * NSTOKES * nsamples );
    data_buffer_vdif  = create_pinned_data_buffer_vdif( nsamples * nchans * npols * npointing * 2 * sizeof(float) );

    if (vm->do_inverse_pfb)
    {
        malloc_ipfb( &gi, filter, nsamples, npols, npointing );
        cu_load_ipfb_filter( filter, &gi );
    }

    // Set up parallel streams
    cudaStream_t streams[npointing];
    for (p = 0; p < npointing; p++)
        cudaStreamCreate(&(streams[p])) ;

    // Create structures for holding header information
    mpi_psrfits mpfs[npointing];
    for (p = 0; p < npointing; p++)
    {
        init_mpi_psrfits(
                &(mpfs[p]),
                vm->obs_metadata,
                vm->vcs_metadata,
                ncoarse_chans,
                mpi_proc_id,
                opts.max_sec_per_file,
                NSTOKES,
                &(beam_geom_vals[p]),
                NULL,
                writer,
                true );
    }

    vdif_header     vhdr;
    struct vdifinfo *vf;
    vf = (struct vdifinfo *)malloc(npointing * sizeof(struct vdifinfo));

    cuDoubleComplex *D = NULL; // See Eqs. (27) to (29) in Ord et al. (2019)

    /****************************
     * GET CALIBRATION SOLUTION *
     ****************************/

    if (cal.cal_type == CAL_RTS)
    {
        D = get_rts_solution( vm->cal_metadata, vm->obs_metadata, cal.caldir, vm->coarse_chan_idxs_to_process[0] );
        if (cal.apply_xy_correction)
        {
            pq_phase_correction( vm->obs_metadata->obs_id, &cal.phase_slope, &cal.phase_offset );

            // Print a suitable message
            if (cal.phase_slope == 0.0 && cal.phase_offset == 0.0)
                logger_timed_message( log, "No XY phase correction information for this obsid" );
            else
            {
                sprintf( log_message, "Applying XY phase correction %.2e*freq%+.2e",
                        cal.phase_slope, cal.phase_offset );
                logger_timed_message( log, log_message );
            }
        }
        else
        {
            logger_timed_message( log, "Not applying XY phase correction" );
        }
    }
    else if (cal.cal_type == CAL_OFFRINGA)
    {
        fprintf( stderr, "error: Offringa-style calibration solutions not currently supported\n" );
        exit(EXIT_FAILURE);
        /*
        // Find the ordering of antennas in Offringa solutions from metafits file
        read_offringa_gains_file( D, nants, cal.offr_chan_num, cal.filename );
        */
    }

    // ------------------
    // Prepare primary beam and geometric delay arrays
    // ------------------
    create_primary_beam( &pb, vm->obs_metadata, vm->coarse_chan_idxs_to_process[0], npointing );
    create_geometric_delays( &gdelays, vm->obs_metadata, vm->vcs_metadata, vm->coarse_chan_idxs_to_process[0], npointing );

    // ------------------

    sprintf( log_message, "Preparing headers for output (receiver channel %lu)",
            vm->obs_metadata->metafits_coarse_chans[vm->coarse_chan_idxs_to_process[0]].rec_chan_number );
    logger_message( log, log_message );

    // Populate the relevant header structs
    populate_vdif_header( vf, &vhdr, vm->obs_metadata, vm->vcs_metadata, vm->coarse_chan_idxs_to_process[0],
            beam_geom_vals, npointing );

    // Begin the main loop: go through data one second at a time

    logger_message( log, "\n*****BEGIN BEAMFORMING*****" );
    char error_message[ERROR_MESSAGE_LEN];

    uintptr_t ntimesteps = vm->num_gps_seconds_to_process;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        gps_second = vm->gps_seconds_to_process[timestep_idx];

        sprintf( log_message, "---Processing GPS second %ld [%lu/%lu]---",
                gps_second, timestep_idx+1, ntimesteps );
        logger_message( log, log_message );

        // Read in data from next file
        logger_start_stopwatch( log, "read", true );

        if (mwalib_voltage_context_read_second(
                    vm->vcs_context,
                    gps_second,
                    1,
                    vm->coarse_chan_idxs_to_process[0],
                    data,
                    data_size,
                    error_message,
                    ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
        {
            fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s", error_message );
            exit(EXIT_FAILURE);
        }
        logger_stop_stopwatch( log, "read" );

        // Get the next second's worth of phases / jones matrices, if needed
        logger_start_stopwatch( log, "delay", true );

        sec_offset = (double)(timestep_idx + vm->gps_seconds_to_process[0] - vm->obs_metadata->obs_id);
        mjd = vm->obs_metadata->sched_start_mjd + (sec_offset + 0.5)/86400.0;
        for (p = 0; p < npointing; p++)
            calc_beam_geom( ras_hours[p], decs_degs[p], mjd, &beam_geom_vals[p] );

        // Calculate the primary beam
        calc_primary_beam( &pb, beam_geom_vals );

        // Calculate the geometric delays
        calc_all_geometric_delays( &gdelays, beam_geom_vals );
        push_geometric_delays_to_device( &gdelays );

        get_jones(
                npointing,              // number of pointings
                vm->obs_metadata,
                vm->coarse_chan_idxs_to_process[0],
                &cal,                   // struct holding info about calibration
                D,                      // Calibration Jones matrices
                pb.B,                   // Primary beam jones matrices
                gf.J );                 // inverse Jones array (output)

        logger_stop_stopwatch( log, "delay" );

        // The writing (of the previous second) is put here in order to
        // allow the possibility that it can overlap with the reading step.
        // Because of this, another "write" has to happen after this loop
        // has terminated
        if (timestep_idx > 0) // i.e. don't do this the first time around
        {
            write_step( mpfs, npointing, vm->output_fine_channels, vm->output_coarse_channels, vf, &vhdr, data_buffer_vdif, log );
        }

        // Form the beams
        logger_start_stopwatch( log, "calc", true );

        cu_form_beam( data, nsamples, gdelays.d_phi, timestep_idx,
                npointing, nants, nchans, npols, invw, &gf,
                detected_beam, data_buffer_coh,
                streams, nchunk, mpfs );

        logger_stop_stopwatch( log, "calc" );

        // Invert the PFB, if requested
        logger_start_stopwatch( log, "ipfb", true );

        if (vm->do_inverse_pfb)
        {
            cu_invert_pfb( detected_beam, timestep_idx, npointing,
                    nsamples, nchans, npols, vf->sizeof_buffer,
                    &gi, data_buffer_vdif );
        }

        logger_stop_stopwatch( log, "ipfb" );

        // Splice channels together
        if (vm->output_fine_channels) // Only PSRFITS output can be combined into a single file
        {
            logger_start_stopwatch( log, "splice", true );

            for (p = 0; p < npointing; p++)
            {
                gather_splice_psrfits( &(mpfs[p]) );
            }

            logger_stop_stopwatch( log, "splice" );
        }
    }

    // Write out the last second's worth of data
    write_step( mpfs, npointing, vm->output_fine_channels, vm->output_coarse_channels, vf, &vhdr, data_buffer_vdif, log );

    logger_message( log, "\n*****END BEAMFORMING*****\n" );

    // Clean up channel-dependent memory
    for (p = 0; p < npointing; p++)
    {
        free_mpi_psrfits( &(mpfs[p]) );

        if (vm->output_coarse_channels)
        {
            free( vf[p].b_scales  );
            free( vf[p].b_offsets );
        }
    }

    logger_report_all_stats( log );
    logger_message( log, "" );

    // Free up memory
    logger_timed_message( log, "Starting clean-up" );

    destroy_detected_beam( detected_beam, npointing, 2*nsamples, nchans );

    free_pfb_filter( filter );

    free( D );

    cudaFreeHost( data_buffer_coh   );
    cudaCheckErrors( "cudaFreeHost(data_buffer_coh) failed" );
    cudaFreeHost( data_buffer_vdif  );
    cudaCheckErrors( "cudaFreeHost(data_buffer_vdif) failed" );
    cudaFreeHost( data );
    cudaCheckErrors( "cudaFreeHost(data) failed" );

    free( opts.pointings_file  );
    free( opts.datadir         );
    free( opts.begin_str       );
    free( opts.coarse_chan_str );
    free( cal.caldir           );
    free( opts.metafits        );
    free( opts.synth_filter    );

    free_formbeam( &gf );
    if (vm->do_inverse_pfb)
    {
        free_ipfb( &gi );
    }

    // Clean up memory associated with the Jones matrices
    free_primary_beam( &pb );
    free_geometric_delays( &gdelays );

    destroy_vcsbeam_metadata( vm );

    // Finalise MPI
    MPI_Finalize();

    return EXIT_SUCCESS;
}


void usage()
{
    printf( "\nusage: make_mwa_tied_array_beam [OPTIONS]\n");

    printf( "\nREQUIRED OPTIONS\n\n"
            "\t-P, --pointings=FILE       FILE containing RA and Decs of multiple pointings\n"
            "\t                           in the format hh:mm:ss.s dd:mm:ss.s ...\n"
            "\t-m, --metafits=FILE        FILE is the metafits file for the target observation\n"
          );

    printf( "\nOPTIONAL OPTIONS\n\n"
            "\t-b, --begin=GPSTIME        Begin time of observation, in GPS seconds\n"
            "\t                           If GPSTIME starts with a '+' or a '-', then the time\n"
            "\t                           is taken relative to the start or end of the observation\n"
            "\t                           respectively. [default: \"+0\"]\n"
            "\t-d, --data-location=PATH   PATH is the directory containing the recombined data\n"
            "\t                           [default: current directory]\n"
            "\t-e, --emulate-legacy       Emulate the legacy system by converting MWAX data into\n"
            "\t                           VCS-recombined-style data using the offline forward fine\n"
            "\t                           PFB. If the input data is already legacy data, this has\n"
            "\t                           no effect [default: off]\n"
            "\t-f, --coarse-chan=CHAN     Coarse channel number\n"
            "\t                           If CHAN starts with a '+' or a '-', then the channel is taken\n"
            "\t                           relative to the first or last channel in the observation\n"
            "\t                           respectively. Otherwise, it is treated as a receiver channel number\n"
            "\t                           (0-255) [default: \"+0\"]\n"
            "\t-p, --out-fine             Output fine-channelised, full-Stokes data (PSRFITS)\n"
            "\t                           (if neither -p nor -v are used, default behaviour is to match channelisation of input)\n"
            "\t-v, --out-coarse           Output coarse-channelised, 2-pol (XY) data (VDIF)\n"
            "\t                           (if neither -p nor -v are used, default behaviour is to match channelisation of input)\n"
            "\t-t, --max_t                Maximum number of seconds per output fits file. [default: 200]\n"
            "\t-S, --synth_filter=filter  Apply the named filter during high-time resolution synthesis.\n"
            "\t                           filter can be MIRROR or LSQ12.\n"
            "\t                           [default: LSQ12]\n"
            "\t-T, --nseconds=VAL         Process VAL seconds of data [default: as many as possible]\n"
          );

    printf( "\nCALIBRATION OPTIONS (RTS)\n\n"
            "\t-c, --cal-metafits=FILE    FILE is the metafits file pertaining to the calibration solution\n"
            "\t-C, --cal-location=PATH    PATH is the directory (RTS) or the file (OFFRINGA) containing the calibration solution\n"
            "\t-R, --ref-ant=ANT          Rotate the phases of the XX and YY elements of the calibration\n"
            "\t                           Jones matrices so that the phases of tile ANT align. If ANT is\n"
            "\t                           outside the range 0-127, no phase rotation is done\n"
            "\t                           [default: 0]\n"
            "\t-X, --cross-terms          Retain the XY and YX terms of the calibration solution\n"
            "\t                           [default: off]\n"
            "\t-U, --no-XY-phase          Do not apply the XY phase correction to the calibration solution\n"
          );

    printf( "\nCALIBRATION OPTIONS (OFFRINGA) -- NOT YET SUPPORTED\n\n"
            "\t-O, --offringa             The calibration solution is in the Offringa format instead of\n"
            "\t                           the default RTS format. In this case, the argument to -C should\n" 
            "\t                           be the full path to the binary solution file.\n"
          );

    printf( "\nOTHER OPTIONS\n\n"
            "\t-h, --help                 Print this help and exit\n"
            "\t-g, --gpu-mem=N            The maximum amount of GPU memory you want make_beam to use in GB [default: -1]\n"
            "\t-V, --version              Print version number and exit\n\n"
          );
}



void make_tied_array_beam_parse_cmdline(
        int argc, char **argv, struct make_tied_array_beam_opts *opts, struct calibration *cal )
{
    // Set defaults
    opts->begin_str          = NULL;  // Absolute or relative GPS time -- when to start beamforming
    opts->nseconds           = -1;    // How many seconds to process (-1 = as many as possible)
    opts->pointings_file     = NULL;  // File containing list of pointings "hh:mm:ss dd:mm:ss ..."
    opts->datadir            = NULL;  // The path to where the recombined data live
    opts->metafits           = NULL;  // filename of the metafits file for the target observation
    opts->coarse_chan_str    = NULL;  // Absolute or relative coarse channel
    opts->out_fine           = false; // Output fine channelised data (PSRFITS)
    opts->out_coarse         = false; // Output coarse channelised data (VDIF)
    opts->emulate_legacy     = false; // Emulate the legacy VCS system using the offline fine PFB
    opts->synth_filter       = NULL;
    opts->max_sec_per_file   = 200;   // Number of seconds per fits files
    opts->gpu_mem            = -1.0;

    cal->metafits            = NULL;  // filename of the metafits file for the calibration observation
    cal->caldir              = NULL;  // The path to where the calibration solutions live
    cal->cal_type            = CAL_RTS;
    cal->ref_ant             = 0;
    cal->cross_terms         = 0;
    cal->phase_offset        = 0.0;
    cal->phase_slope         = 0.0;
    cal->apply_xy_correction = true;

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"out-fine",        no_argument,       0, 'p'},
                {"out-coarse",      no_argument,       0, 'v'},
                {"emulate_legacy",  no_argument,       0, 'e'},
                {"max_t",           required_argument, 0, 't'},
                {"synth_filter",    required_argument, 0, 'S'},
                {"nseconds",        required_argument, 0, 'T'},
                {"pointings",       required_argument, 0, 'P'},
                {"data-location",   required_argument, 0, 'd'},
                {"cal-location",    required_argument, 0, 'C'},
                {"metafits",        required_argument, 0, 'm'},
                {"cal-metafits",    required_argument, 0, 'c'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"ref-ant",         required_argument, 0, 'R'},
                {"cross-terms",     no_argument,       0, 'X'},
                {"no-XY-phase",     required_argument, 0, 'U'},
                {"offringa",        no_argument      , 0, 'O'},
                {"gpu-mem",         required_argument, 0, 'g'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "b:c:C:d:e:f:F:g:hm:OpP:R:S:t:T:UvVX",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c)
            {
                case 'b':
                    opts->begin_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->begin_str, optarg );
                    break;
                case 'c':
                    cal->metafits = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( cal->metafits, optarg );
                    break;
                case 'C':
                    cal->caldir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( cal->caldir, optarg );
                    break;
                case 'd':
                    opts->datadir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->datadir, optarg );
                    break;
                case 'e':
                    opts->emulate_legacy = true;
                    break;
                case 'f':
                    opts->coarse_chan_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->coarse_chan_str, optarg );
                    break;
                case 'g':
                    opts->gpu_mem = atof(optarg);
                    break;
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'O':
                    cal->cal_type = CAL_OFFRINGA;
                    break;
                case 'p':
                    opts->out_fine = true;
                    break;
                case 'P':
                    opts->pointings_file = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->pointings_file, optarg );
                    break;
                case 'R':
                    cal->ref_ant = atoi(optarg);
                    break;
                case 'S':
                    opts->synth_filter = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->synth_filter, optarg );
                    break;
                case 't':
                    opts->max_sec_per_file = atoi(optarg);
                    break;
                case 'T':
                    opts->nseconds = atol(optarg);
                    if (opts->nseconds <= 0)
                    {
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "-%c argument must be >= 1\n", c );
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'U':
                    cal->apply_xy_correction = false;
                    break;
                case 'v':
                    opts->out_coarse = true;
                    break;
                case 'V':
                    printf( "MWA Beamformer %s\n", VCSBEAM_VERSION);
                    exit(0);
                    break;
                case 'X':
                    cal->cross_terms = 1;
                    break;
                default:
                    fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                    "unrecognised option '%s'\n", optarg );
                    usage();
                    exit(EXIT_FAILURE);
            }
        }
    }
    else
    {
        usage();
        exit(EXIT_FAILURE);
    }


    // Check that all the required options were supplied
    assert( opts->pointings_file != NULL );
    assert( cal->caldir          != NULL );
    assert( opts->metafits       != NULL );
    assert( cal->metafits        != NULL );

    if (opts->datadir == NULL)
    {
        opts->datadir = (char *)malloc( 2 );
        strcpy( opts->datadir, "." );
    }

    if (opts->begin_str == NULL)
    {
        opts->begin_str = (char *)malloc( 3 );
        strcpy( opts->begin_str, "+0" );
    }

    if (opts->coarse_chan_str == NULL)
    {
        opts->coarse_chan_str = (char *)malloc( 3 );
        strcpy( opts->coarse_chan_str, "+0" );
    }
}


void parse_pointing_file( const char *filename, double **ras_hours, double **decs_degs, unsigned int *npointings )
/* Parse the given file in FILENAME and create arrays of RAs and DECs.
 * This function allocates memory for ras_hours and decs_degs arrays.
 * Caller can destroy with free().
 */
{
    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error: cannot open pointings file %s\n", filename );
        exit(EXIT_FAILURE);
    }

    // Do one pass through the file to count "words"
    // The RAs and Decs are expected to be whitespace-delimited
    int nwords = 0;
    char word[64];
    while (fscanf( f, "%s", word ) != EOF)
        nwords++;

    // Check that we have an even number of words (they should be in RA/Dec pairs)
    if (nwords % 2 != 0)
    {
        fprintf( stderr, "error: cannot parse pointings file %s\n", filename );
        exit(EXIT_FAILURE);
    }
    *npointings = nwords/2;

    // Allocate memory
    *ras_hours = (double *)malloc( *npointings * sizeof(double) );
    *decs_degs = (double *)malloc( *npointings * sizeof(double) );

    // Rewind to beginning of file, read the words in again, and parse them
    rewind( f );
    char ra_str[64], dec_str[64];
    unsigned int p;
    for (p = 0; p < *npointings; p++)
    {
        // Read in the next Ra/Dec pair
        fscanf( f, "%s %s", ra_str, dec_str );

        // Parse them and make them decimal
        (*ras_hours)[p] = parse_ra( ra_str );
        (*decs_degs)[p] = parse_dec( dec_str );
    }

    // Close the file
    fclose( f );
}


void write_step( mpi_psrfits *mpfs, int npointing, bool out_fine, bool out_coarse,
        struct vdifinfo *vf, vdif_header *vhdr, float *data_buffer_vdif, logger *log )
{
    int mpi_proc_id;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_proc_id );

    int p;
    for (p = 0; p < npointing; p++)
    {
        logger_start_stopwatch( log, "write", true );

        if (out_fine)
        {
            if (mpi_proc_id == mpfs[p].writer_id)
            {
                // Write out the spliced channels
                wait_splice_psrfits( &(mpfs[p]) );

                if (psrfits_write_subint( &(mpfs[p].spliced_pf) ) != 0)
                {
                    fprintf(stderr, "error: Write PSRFITS subint failed. File exists?\n");
                    exit(EXIT_FAILURE);
                }
                mpfs[p].spliced_pf.sub.offs = roundf(mpfs[p].spliced_pf.tot_rows * mpfs[p].spliced_pf.sub.tsubint) + 0.5*mpfs[p].spliced_pf.sub.tsubint;
                mpfs[p].spliced_pf.sub.lst += mpfs[p].spliced_pf.sub.tsubint;
            }
        }

        if (out_coarse)
        {
            vdif_write_second( &vf[p], vhdr,
                    data_buffer_vdif + p * vf->sizeof_buffer );
        }

        logger_stop_stopwatch( log, "write" );
    }
}
