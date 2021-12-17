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
    char              *custom_flags;    // Use custom list for flagging antennas

    // Output options
    bool               out_fine;         // Output fine channelised data (PSRFITS)
    bool               out_coarse;       // Output coarse channelised data (VDIF)
    bool               emulate_legacy;   // Use the offline forward pfb if input data is MWAX

    // Other options
    char              *synth_filter;     // Which synthesis filter to use
    int                max_sec_per_file; // Number of seconds per fits files
    float              gpu_mem_GB;       // Default = 0.0. If 0.0 use all GPU mem
};

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/

void usage();
void make_tied_array_beam_parse_cmdline( int argc, char **argv, struct make_tied_array_beam_opts *opts, calibration *cal );

void write_step( vcsbeam_context *vm, mpi_psrfits *mpfs,
        struct vdifinfo *vf, vdif_header *vhdr, float *data_buffer_vdif );

/********
 * MAIN *
 ********/

int main(int argc, char **argv)
{
    // Parse command line arguments
    struct make_tied_array_beam_opts opts;
    calibration cal;           // Variables for calibration settings
    init_calibration( &cal );
    make_tied_array_beam_parse_cmdline( argc, argv, &opts, &cal );

    // Initialise MPI
    MPI_Init( NULL, NULL );
    int world_size, mpi_proc_id;
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_proc_id );

    const int writer = 0; // Designate process 0 to write out the files
    const int ncoarse_chans = world_size;

    int i; // Generic counter

    // Create an mwalib metafits context and associated metadata
    const int chans_per_proc = 1;
    vcsbeam_context *vm = init_vcsbeam_context(
        opts.metafits, cal.metafits,
        opts.coarse_chan_str, chans_per_proc, mpi_proc_id,
        opts.begin_str, opts.nseconds, 0,
        opts.datadir );

    char mwalib_version[32];
    get_mwalib_version( mwalib_version );

    sprintf( vm->log_message, "Creating metafits and voltage contexts via MWALIB (v%s)",
            mwalib_version );
    logger_timed_message( vm->log, vm->log_message );

    if (opts.out_fine)    set_vcsbeam_fine_output( vm, true );
    if (opts.out_coarse)  set_vcsbeam_coarse_output( vm, true );

    // Create some "shorthand" variables for code brevity
    uintptr_t nchans         = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t npols          = vm->obs_metadata->num_ant_pols;   // (X,Y)
    unsigned int nsamples    = vm->sample_rate;

    // Parse input pointings
    vmParsePointingFile( vm, opts.pointings_file );

    // Allocate memory for various data products
    cuDoubleComplex  ****detected_beam = create_detected_beam( vm->npointing, 2*nsamples, nchans, npols );

    // Get pointing geometry information
    beam_geom beam_geom_vals[vm->npointing];

    unsigned int p;
    double mjd, sec_offset;
    mjd = vm->obs_metadata->sched_start_mjd;
    for (p = 0; p < vm->npointing; p++)
        calc_beam_geom( vm->ras_hours[p], vm->decs_degs[p], mjd, &beam_geom_vals[p] );

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

    /* Allocate host and device memory for the use of the cu_form_beam function */
    // Declaring pointers to the structs so the memory can be alternated
    vmSetMaxGPUMem( vm, (uintptr_t)(opts.gpu_mem_GB * 1024*1024*1024) );
    vmMallocVHost( vm );
    vmMallocVDevice( vm );
    vmMallocJVDevice( vm );
    vmMallocEHost( vm );
    vmMallocEDevice( vm );
    vmMallocSHost( vm );
    vmMallocSDevice( vm );
    vmMallocJHost( vm );
    vmMallocJDevice( vm );
    vmMallocDHost( vm );
    vmMallocDDevice( vm );
    vmMallocPQIdxsHost( vm );
    vmMallocPQIdxsDevice( vm );

    // Create a lists of rf_input indexes ordered by antenna number (needed for gpu kernels)
    vmSetPolIdxLists( vm );

    // ... and upload them to the gpu, ready for use!
    vmPushPolIdxLists( vm );

    // Create output buffer arrays

    struct gpu_ipfb_arrays gi;
    float *data_buffer_vdif   = NULL;
    if (vm->do_inverse_pfb)
    {
        data_buffer_vdif  = create_pinned_data_buffer( nsamples * nchans * npols * vm->npointing * 2 * sizeof(float) );
        malloc_ipfb( &gi, filter, nsamples, npols, vm->npointing );
        cu_load_ipfb_filter( filter, &gi );
    }

    // Create structures for holding header information
    mpi_psrfits mpfs[vm->npointing];
    for (p = 0; p < vm->npointing; p++)
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
    vf = (struct vdifinfo *)malloc(vm->npointing * sizeof(struct vdifinfo));

    /****************************
     * GET CALIBRATION SOLUTION *
     ****************************/

    if (cal.cal_type == CAL_RTS)
    {
        vmLoadRTSSolution( vm, cal.use_bandpass, cal.caldir, mpi_proc_id );
    }
    else if (cal.cal_type == CAL_OFFRINGA)
    {
        vmLoadOffringaSolution( vm, mpi_proc_id, cal.caldir );
    }

    // Flag antennas that need flagging
    vmSetCustomTileFlags( vm, opts.custom_flags, &cal );

    // Apply any calibration corrections
    parse_calibration_correction_file( vm->obs_metadata->obs_id, &cal );
    apply_calibration_corrections( &cal, vm->D, vm->obs_metadata,
            vm->coarse_chan_idxs_to_process[0], vm->log );

    // ------------------
    // Prepare primary beam and geometric delay arrays
    // ------------------
    vmCreatePrimaryBeam( vm );
    vmCreateGeometricDelays( vm );
    vmCreateStatistics( vm, mpfs );

    // ------------------

    sprintf( vm->log_message, "Preparing headers for output (receiver channel %lu)",
            vm->obs_metadata->metafits_coarse_chans[vm->coarse_chan_idxs_to_process[0]].rec_chan_number );
    logger_message( vm->log, vm->log_message );

    // Populate the relevant header structs
    populate_vdif_header( vf, &vhdr, vm->obs_metadata, vm->vcs_metadata, vm->coarse_chan_idxs_to_process[0],
            beam_geom_vals, vm->npointing );

    // Begin the main loop: go through data one second at a time

    logger_message( vm->log, "\n*****BEGIN BEAMFORMING*****" );

    uintptr_t ntimesteps = vm->num_gps_seconds_to_process;
    uintptr_t timestep_idx;

    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        // Read in data from next file
        vmReadNextSecond( vm );

        // Calculate J (inverse) and Phi (geometric delays)
        vmCalcJonesAndDelays( vm, vm->ras_hours, vm->decs_degs, beam_geom_vals );

        // Move the needed (just calculated) quantities to the GPU
        vmPushPhi( vm );
        vmPushJ( vm );

        // The writing (of the previous second) is put here in order to
        // allow the possibility that it can overlap with the reading step.
        // Because of this, another "write" has to happen after this loop
        // has terminated
        if (timestep_idx > 0) // i.e. don't do this the first time around
        {
            write_step( vm, mpfs, vf, &vhdr, data_buffer_vdif );
        }

        // Form the beams
        logger_start_stopwatch( vm->log, "calc", true );

        vmBeamformSecond( vm );
        vmPullS( vm );
        vmSendSToFits( vm, mpfs );

        logger_stop_stopwatch( vm->log, "calc" );

        // Invert the PFB, if requested
        if (vm->do_inverse_pfb)
        {
            logger_start_stopwatch( vm->log, "ipfb", true );

            vmPullE( vm );
            prepare_detected_beam( detected_beam, mpfs, vm );
            cu_invert_pfb( detected_beam, timestep_idx, vm->npointing,
                    nsamples, nchans, npols, vf->sizeof_buffer,
                    &gi, data_buffer_vdif );

            logger_stop_stopwatch( vm->log, "ipfb" );
        }

        // Splice channels together
        if (vm->output_fine_channels) // Only PSRFITS output can be combined into a single file
        {
            logger_start_stopwatch( vm->log, "splice", true );

            for (p = 0; p < vm->npointing; p++)
            {
                gather_splice_psrfits( &(mpfs[p]) );
            }

            logger_stop_stopwatch( vm->log, "splice" );
        }
    }

    // Write out the last second's worth of data
    write_step( vm, mpfs, vf, &vhdr, data_buffer_vdif );

    logger_message( vm->log, "\n*****END BEAMFORMING*****\n" );

    // Clean up channel-dependent memory
    for (p = 0; p < vm->npointing; p++)
    {
        free_mpi_psrfits( &(mpfs[p]) );

        if (vm->output_coarse_channels)
        {
            free( vf[p].b_scales  );
            free( vf[p].b_offsets );
        }
    }

    logger_report_all_stats( vm->log );
    logger_message( vm->log, "" );

    // Free up memory
    logger_timed_message( vm->log, "Starting clean-up" );

    destroy_detected_beam( detected_beam, vm->npointing, 2*nsamples, nchans );

    free_pfb_filter( filter );

    if (vm->do_inverse_pfb)
    {
        cudaFreeHost( data_buffer_vdif  );
        cudaCheckErrors( "cudaFreeHost(data_buffer_vdif) failed" );
    }

    vmDestroyStatistics( vm );

    free( opts.pointings_file  );
    free( opts.datadir         );
    free( opts.begin_str       );
    free( opts.coarse_chan_str );
    free( opts.custom_flags    );
    free( opts.metafits        );
    free( opts.synth_filter    );

    free_calibration( &cal );

    vmFreeVHost( vm );
    vmFreeVDevice( vm );
    vmFreeJVDevice( vm );
    vmFreeEHost( vm );
    vmFreeEDevice( vm );
    vmFreeSHost( vm );
    vmFreeSDevice( vm );
    vmFreeJDevice( vm );
    vmFreeJHost( vm );
    vmFreeDHost( vm );
    vmFreeDDevice( vm );
    vmFreePQIdxsDevice( vm );
    vmFreePQIdxsHost( vm );


    if (vm->do_inverse_pfb)
    {
        free_ipfb( &gi );
    }

    // Clean up memory associated with the Jones matrices
    free_primary_beam( &vm->pb );
    free_geometric_delays( &vm->gdelays );

    // Free the CUDA streams
    vmDestroyCudaStreams( vm );

    // Destroy the whole context
    destroy_vcsbeam_context( vm );

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
            "\t-B, --bandpass             Use the Bandpass calibrations (as well as the DIJones solutions) [default: off]\n"
            "\t-c, --cal-metafits=FILE    FILE is the metafits file pertaining to the calibration solution\n"
            "\t-C, --cal-location=PATH    PATH is the directory (RTS) or the file (OFFRINGA) containing the calibration solution\n"
            "\t-F  --flagged-tiles=FILE   FILE is a text file including the TileNames of extra tiles to be flagged.\n"
            "\t                           By default, tiles flagged in both the calibration and the observation metafits file\n"
            "\t                           are flagged in the beamformer. The utility 'rts_flag_ant_to_tilenames.py' can be used\n"
            "\t                           to convert the antenna numbers listed in the RTS 'flagged_tiles.txt' file into human-\n"
            "\t                           readable tile names.\n"
            "\t-R, --ref-ant=TILENAME     Override the reference tile given in pq_phase_correction.txt for rotating the phases\n"
            "\t                           of the PP and QQ elements of the calibration solution. To turn off phase rotation\n"
            "\t                           altogether, set TILENAME=NONE.\n"
            "\t-U, --PQ-phase=PH,OFFS     Override the phase correction given in pq_phase_correction.txt. PH is given in rad/Hz\n"
            "\t                           and OFFS given in rad, such that, the QQ element of the calibration Jones matrix\n"
            "\t                           for frequency F (in Hz) is multiplied by\n"
            "\t                                exp(PH*F + OFFS)\n"
            "\t                           Setting PH = OFFS = 0 is equivalent to not performing any phase correction\n"
            "\t-X, --cross-terms          Retain the PQ and QP terms of the calibration solution [default: off]\n"
          );

    printf( "\nCALIBRATION OPTIONS (OFFRINGA)\n\n"
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
        int argc, char **argv, struct make_tied_array_beam_opts *opts, calibration *cal )
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
    opts->gpu_mem_GB         = 0.0;
    opts->custom_flags       = NULL;

    cal->metafits            = NULL;  // filename of the metafits file for the calibration observation
    cal->caldir              = NULL;  // The path to where the calibration solutions live
    cal->cal_type            = CAL_RTS;
    cal->ref_ant             = NULL;
    cal->keep_cross_terms    = false;
    cal->phase_offset        = 0.0;
    cal->phase_slope         = 0.0;
    cal->custom_pq_correction = false;
    cal->use_bandpass        = false; // use the Bandpass calibration solutions

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"bandpass",        no_argument,       0, 'B'},
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
                {"flagged-tiles",   required_argument, 0, 'F'},
                {"cross-terms",     no_argument,       0, 'X'},
                {"PQ-phase",        required_argument, 0, 'U'},
                {"offringa",        no_argument      , 0, 'O'},
                {"gpu-mem",         required_argument, 0, 'g'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "b:Bc:C:d:e:f:F:g:hm:OpP:R:S:t:T:U:vVX",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c)
            {
                case 'b':
                    opts->begin_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->begin_str, optarg );
                    break;
                case 'B':
                    cal->use_bandpass = true;
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
                case 'F':
                    opts->custom_flags = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->custom_flags, optarg );
                    break;
                case 'g':
                    opts->gpu_mem_GB = atof(optarg);
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
                    cal->ref_ant = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( cal->ref_ant, optarg );
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
                    if (sscanf( optarg, "%lf,%lf", &(cal->phase_slope), &(cal->phase_offset) ) != 2)
                    {
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "cannot parse -U option (\"%s\") as \"FLOAT,FLOAT\"\n", optarg );
                        exit(EXIT_FAILURE);
                    }
                    cal->custom_pq_correction = true;
                    break;
                case 'v':
                    opts->out_coarse = true;
                    break;
                case 'V':
                    printf( "MWA Beamformer %s\n", VCSBEAM_VERSION);
                    exit(0);
                    break;
                case 'X':
                    cal->keep_cross_terms = true;
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



void write_step( vcsbeam_context *vm, mpi_psrfits *mpfs,
        struct vdifinfo *vf, vdif_header *vhdr, float *data_buffer_vdif )
{
    int mpi_proc_id;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_proc_id );

    int p;
    for (p = 0; p < vm->npointing; p++)
    {
        logger_start_stopwatch( vm->log, "write", true );

        if (vm->output_fine_channels)
        {
            // Write out the spliced channels
            wait_splice_psrfits( &(mpfs[p]) );

            if (mpi_proc_id == mpfs[p].writer_id)
            {
                if (psrfits_write_subint( &(mpfs[p].spliced_pf) ) != 0)
                {
                    fprintf(stderr, "error: Write PSRFITS subint failed. File exists?\n");
                    exit(EXIT_FAILURE);
                }
                mpfs[p].spliced_pf.sub.offs = roundf(mpfs[p].spliced_pf.tot_rows * mpfs[p].spliced_pf.sub.tsubint) + 0.5*mpfs[p].spliced_pf.sub.tsubint;
                mpfs[p].spliced_pf.sub.lst += mpfs[p].spliced_pf.sub.tsubint;
            }
        }

        if (vm->output_coarse_channels)
        {
            vdif_write_second( &vf[p], vhdr,
                    data_buffer_vdif + p * vf->sizeof_buffer );
        }

        logger_stop_stopwatch( vm->log, "write" );
    }
}
