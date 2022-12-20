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
    char              *custom_flags;     // Use custom list for flagging antennas

    // Output format options
    bool               out_psrfits;      // Output PSRFITS format
    bool               out_vdif;         // Output VDIF format
    bool               out_filterbank;   // Output FILTERBANK format

    // These are worked out from the output format options
    bool               out_fine;         // Output fine channelised data
    bool               out_coarse;       // Output coarse channelised data

    // Calibration options
    char              *cal_metafits;     // Filename of the metafits file
    char              *caldir;           // Location of calibration data
    int                cal_type;         // Either RTS or OFFRINGA
    char              *ref_ant;          // Reference antenna for calibration phases
    double             phase_offset;     // Rotate the phase of Y by m*freq + c, where
    double             phase_slope;      //   m = phase_slope (rad/Hz)
                                         //   c = phase_offset (rad)
    bool               custom_pq_correction; // Set to true if phase_offset and phase_slope are to be applied
    bool               keep_cross_terms; // Include PQ and QP of calibration Jones matrices
    bool               use_bandpass;     // Use the Bandpass solutions

    // Other options
    char              *analysis_filter;  // Which analysis filter to use
    char              *synth_filter;     // Which synthesis filter to use
    bool               smart;            // Use legacy settings for PFB
    int                max_sec_per_file; // Number of seconds per fits files
    //float              gpu_mem_GB;       // Default = 0.0. If 0.0 use all GPU mem
    int                nchunks;          // Split each second into this many processing chunks
};

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/

void usage();
void make_tied_array_beam_parse_cmdline( int argc, char **argv, struct make_tied_array_beam_opts *opts );

void write_step( vcsbeam_context *vm,
#ifdef HAVE_PSRFITS
        mpi_psrfits *mpfs,
#endif
#ifdef HAVE_VDIF
        struct vdifinfo *vf, vdif_header *vhdr, float *data_buffer_vdif,
#endif
        struct make_tied_array_beam_opts *opts );

/********
 * MAIN *
 ********/

int main(int argc, char **argv)
{
    // Parse command line arguments
    struct make_tied_array_beam_opts opts;
    make_tied_array_beam_parse_cmdline( argc, argv, &opts );

    int i; // Generic counter

    // Create an mwalib metafits context and associated metadata
    bool use_mpi = true;
    vcsbeam_context *vm = vmInit( use_mpi );

    vmPrintTitle( vm, "Beamformer" );

    vmLoadObsMetafits( vm, opts.metafits );
    vmLoadCalMetafits( vm, opts.cal_metafits );

    vmBindObsData( vm,
        opts.coarse_chan_str, 1, vm->mpi_rank,
        opts.begin_str, opts.nseconds, 0,
        opts.datadir );

    // If explicit output flags are given, set the output channelisation
    // accordingly
    if (opts.out_fine || opts.out_coarse)
    {
        vmSetOutputChannelisation( vm, opts.out_fine, opts.out_coarse );
    }

    // If we need to, set up the forward PFB
    vm->chunks_per_second = opts.nchunks;
    if (vm->do_forward_pfb)
    {
        // Load the filter
        int K = 128; // The number of desired output channels
        vmLoadFilter( vm, opts.analysis_filter, ANALYSIS_FILTER, K );

        // Create and init the PFB struct
        int M = K; // The filter stride (M = K <=> "critically sampled PFB")
        vmInitForwardPFB( vm, M, (opts.smart ? PFB_SMART : PFB_FULL_PRECISION) );
    }

    vm->cal.metafits     = strdup( opts.cal_metafits );
    vm->cal.ref_ant      = strdup( opts.ref_ant );
    vm->cal.phase_offset = opts.phase_offset;
    vm->cal.phase_slope  = opts.phase_slope;
    vm->cal.custom_pq_correction = opts.custom_pq_correction;
    vm->cal.keep_cross_terms     = opts.keep_cross_terms;

    char mwalib_version[32];
    get_mwalib_version( mwalib_version );

    sprintf( vm->log_message, "Creating metafits and voltage contexts via MWALIB (v%s)",
            mwalib_version );
    logger_timed_message( vm->log, vm->log_message );

    // Parse input pointings
    vmParsePointingFile( vm, opts.pointings_file );

    // Get pointing geometry information
    beam_geom beam_geom_vals[vm->npointing];

    unsigned int p;
    double mjd, sec_offset;
    mjd = vm->obs_metadata->sched_start_mjd;
    for (p = 0; p < vm->npointing; p++)
        calc_beam_geom( vm->ras_hours[p], vm->decs_degs[p], mjd, &beam_geom_vals[p] );

    // Some shorthand variables (only needed for things relating to the inverse PFB)
    uintptr_t nchans         = vm->nfine_chan;
    uintptr_t npols          = vm->obs_metadata->num_ant_pols;
    unsigned int nsamples    = vm->fine_sample_rate;

    cuDoubleComplex  ****detected_beam;

    if (vm->do_inverse_pfb)
    {
        // Load the (synthesis) PFB filter:
        // Adjust by the scaling that was introduced by the forward PFB,
        // along with any other scaling that I, Lord and Master of the inverse
        // PFB, feel is appropriate.
        vmLoadFilter( vm, opts.synth_filter, SYNTHESIS_FILTER, nchans );
        vmScaleFilterCoeffs( vm, SYNTHESIS_FILTER, 15.0/7.2 ); // (1/7.2) = 16384/117964.8

        // Allocate memory for various data products
        detected_beam = create_detected_beam( vm->npointing, 2*nsamples, nchans, npols );
    }

    /*********************
     * Memory Allocation *
     *********************/

    /* Allocate host and device memory for the use of the cu_form_beam function */
    // Declaring pointers to the structs so the memory can be alternated
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
    // and upload them to the gpu
    vmSetPolIdxLists( vm );
    vmPushPolIdxLists( vm );

    // Create output buffer arrays

    struct gpu_ipfb_arrays gi;
#ifdef HAVE_VDIF
    float *data_buffer_vdif   = NULL;
    if (vm->do_inverse_pfb)
    {
        data_buffer_vdif  = create_pinned_data_buffer( nsamples * nchans * npols * vm->npointing * 2 * sizeof(float) );
        malloc_ipfb( &gi, vm->synth_filter, nsamples, npols, vm->npointing );
        cu_load_ipfb_filter( vm->synth_filter, &gi );
    }
#endif

    // Create structures for holding header information
#ifdef HAVE_PSRFITS
    mpi_psrfits mpfs[vm->npointing];
    if (opts.out_psrfits)
    {
        for (p = 0; p < vm->npointing; p++)
        {
            vmInitMPIPsrfits( vm, &(mpfs[p]), opts.max_sec_per_file, NSTOKES,
                    &(beam_geom_vals[p]), NULL, true );
        }
    }
#endif

    /****************************
     * GET CALIBRATION SOLUTION *
     ****************************/

    vmBindCalibrationData( vm, opts.caldir, opts.cal_type, opts.use_bandpass, opts.custom_flags );
    vmReadCalibration( vm );

    // Apply any calibration corrections
    parse_calibration_correction_file( vm->obs_metadata->obs_id, &vm->cal );
    vmApplyCalibrationCorrections( vm );

    // ------------------
    // Prepare primary beam and geometric delay arrays
    // ------------------
    vmCreatePrimaryBeam( vm );
    vmCreateGeometricDelays( vm );
#ifdef HAVE_PSRFITS
    vmCreateStatistics( vm, mpfs );
#endif

    // ------------------

#ifdef HAVE_VDIF
    // Populate the relevant header structs
    vmPopulateVDIFHeader( vm, beam_geom_vals );
#endif

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
            write_step(
                    vm,
#ifdef HAVE_PSRFITS
                    mpfs,
#endif
#ifdef HAVE_VDIF
                    vm->vf, &vm->vhdr, data_buffer_vdif,
#endif
                    &opts );
        }

        // Do the forward PFB (if needed), and form the beams
        vmBeamformSecond( vm );

        // Invert the PFB, if requested
        if (vm->do_inverse_pfb)
        {
            logger_start_stopwatch( vm->log, "ipfb", true );

            vmPullE( vm );
            prepare_detected_beam( detected_beam, vm );
#ifdef HAVE_VDIF
// The need for this HAVE_VDIF counts as a bug. The inverse PFB should work
// without knowing anything about the VDIF output. However, cu_invert_pfb()
// currently uses the sizeof_buffer from the vdif struct.
            cu_invert_pfb( detected_beam, timestep_idx, vm->npointing,
                    nsamples, nchans, npols, vm->vf->sizeof_buffer,
                    &gi, data_buffer_vdif );
#endif

            logger_stop_stopwatch( vm->log, "ipfb" );
        }

#ifdef HAVE_PSRFITS
        // Splice channels together for PSRFITS
        if (opts.out_psrfits)
        {
            vmPullS( vm );
            vmSendSToFits( vm, mpfs );

            logger_start_stopwatch( vm->log, "splice", true );

            for (p = 0; p < vm->npointing; p++)
            {
                gather_splice_psrfits( &(mpfs[p]) );
            }

            logger_stop_stopwatch( vm->log, "splice" );
        }
#endif
    }

    // Write out the last second's worth of data
    write_step(
            vm,
#ifdef HAVE_PSRFITS
            mpfs,
#endif
#ifdef HAVE_VDIF
            vm->vf, &vm->vhdr, data_buffer_vdif,
#endif
            &opts );

    logger_message( vm->log, "\n*****END BEAMFORMING*****\n" );

    // Clean up channel-dependent memory
    for (p = 0; p < vm->npointing; p++)
    {
#ifdef HAVE_PSRFITS
        if (opts.out_psrfits)
        {
            free_mpi_psrfits( &(mpfs[p]) );
        }
#endif

#ifdef HAVE_VDIF
        if (opts.out_vdif)
        {
            free( vm->vf[p].b_scales  );
            free( vm->vf[p].b_offsets );
        }
#endif
    }

    // Report performace statistics
    vmReportPerformanceStats( vm );

    // Free up memory
    logger_timed_message( vm->log, "Starting clean-up" );

#ifdef HAVE_VDIF
    if (vm->do_inverse_pfb)
    {
        destroy_detected_beam( detected_beam, vm->npointing, 2*nsamples, nchans );
    }

    if (vm->do_inverse_pfb)
    {
        cudaFreeHost( data_buffer_vdif  );
        cudaCheckErrors( "cudaFreeHost(data_buffer_vdif) failed" );
    }
#endif

#ifdef HAVE_PSRFITS
    vmDestroyStatistics( vm );
#endif

    free( opts.pointings_file  );
    free( opts.datadir         );
    free( opts.begin_str       );
    free( opts.coarse_chan_str );
    free( opts.custom_flags    );
    free( opts.metafits        );
    free( opts.synth_filter    );

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

    return EXIT_SUCCESS;
}


void usage()
{
    printf( "\nusage: make_mwa_tied_array_beam [OPTIONS]\n");

    printf( "\nREQUIRED OPTIONS\n\n"
            "\t-c, --cal-metafits=FILE    FILE is the metafits file pertaining to the calibration solution\n"
            "\t-C, --cal-location=PATH    PATH is the directory (RTS) or the file (OFFRINGA) containing the calibration solution\n"
            "\t-m, --metafits=FILE        FILE is the metafits file for the target observation\n"
            "\t-P, --pointings=FILE       FILE containing RA and Decs of multiple pointings\n"
            "\t                           in the format hh:mm:ss.s dd:mm:ss.s ...\n"
          );

    printf( "\nCHANNELISATION OPTIONS\n\n"
            "\t-A, --analysis_filter=FILTER  Apply the named filter during fine channelisation (for MWAX only).\n"
            "\t                           [default: FINEPFB]\n"
            "\t-s, --smart                Use legacy settings for fine channelisation [default: off]\n"
            "\t-S, --synth_filter=FILTER  Apply the named filter during high-time resolution synthesis.\n"
            "\t                           FILTER can be MIRROR or LSQ12.\n"
            "\t                           [default: LSQ12]\n"
          );

    printf( "\nINPUT OPTIONS\n\n"
            "\t-b, --begin=GPSTIME        Begin time of observation, in GPS seconds\n"
            "\t                           If GPSTIME starts with a '+' or a '-', then the time\n"
            "\t                           is taken relative to the start or end of the observation\n"
            "\t                           respectively. [default: \"+0\"]\n"
            "\t-d, --data-location=PATH   PATH is the directory containing the recombined data\n"
            "\t                           [default: current directory]\n"
            "\t-f, --coarse-chan=CHAN     Coarse channel number\n"
            "\t                           If CHAN starts with a '+' or a '-', then the channel is taken\n"
            "\t                           relative to the first or last channel in the observation\n"
            "\t                           respectively. Otherwise, it is treated as a receiver channel number\n"
            "\t                           (0-255) [default: \"+0\"]\n"
            "\t-T, --nseconds=VAL         Process VAL seconds of data [default: as many as possible]\n"
          );

    printf( "\nOUTPUT OPTIONS\n\n"
            "\t-o, --output-format=FMT    Define the output format. To generate more than one output format, \n"
            "\t                           repeat -o for each format. At least one output format must be supplied.\n"
            "\t                           Available formats (FMT):\n"
#ifdef HAVE_PSRFITS
            "\t                               PSRFITS    (  10 kHz channels) ✓\n"
#else
            "\t                               PSRFITS    (  10 kHz channels) ✗ (requires PSRFITS_UTILS)\n"
#endif
#ifdef HAVE_VDIF
            "\t                               FILTERBANK (  10 kHz channels) ✓\n"
#else
            "\t                               FILTERBANK (  10 kHz channels) ✗ (requires SIGPROC)\n"
#endif
#ifdef HAVE_VDIF
            "\t                               VDIF       (1.28 MHz channels) ✓\n"
#else
            "\t                               VDIF       (1.28 MHz channels) ✗ (requires VDIFIO)\n"
#endif
            "\t-t, --max_t                Maximum number of seconds per output PSRFITS file. [default: 200]\n"
          );

    printf( "\nCALIBRATION OPTIONS\n\n"
            "\t-B, --bandpass             Use the Bandpass (fine channel) as well as the DIJones (coarse channel) solutions\n"
            "\t                           (only relevant for RTS) [default: off]\n"
            "\t-F, --flagged-tiles=FILE   FILE is a text file including the TileNames of extra tiles to be flagged.\n"
            "\t                           By default, tiles flagged in both the calibration and the observation metafits file\n"
            "\t                           are flagged in the beamformer. The utility 'rts_flag_ant_to_tilenames.py' can be used\n"
            "\t                           to convert the antenna numbers listed in the RTS 'flagged_tiles.txt' file into human-\n"
            "\t                           readable tile names.\n"
            "\t-O, --offringa             The calibration solution is in the Offringa format instead of\n"
            "\t                           the default RTS format. In this case, the argument to -C should\n" 
            "\t                           be the full path to the binary solution file.\n"
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

    printf( "\nMEMORY OPTIONS\n\n"
            "\t-n, --nchunks=VAL          Split each second's worth of data into VAL processing chunks\n"
            "\t                           [default: 1]\n"
          );

    printf( "\nOTHER OPTIONS\n\n"
            "\t-h, --help                 Print this help and exit\n"
            "\t-V, --version              Print version number and exit\n\n"
          );
// This option is currently too problematic to deal with:
//            "\t-g, --gpu-mem=N            The maximum amount of GPU memory you want make_beam to use in GB [default: -1]\n"
}



void make_tied_array_beam_parse_cmdline(
        int argc, char **argv, struct make_tied_array_beam_opts *opts )
{
    // Set defaults
    opts->begin_str            = NULL;  // Absolute or relative GPS time -- when to start beamforming
    opts->nseconds             = -1;    // How many seconds to process (-1 = as many as possible)
    opts->pointings_file       = NULL;  // File containing list of pointings "hh:mm:ss dd:mm:ss ..."
    opts->datadir              = NULL;  // The path to where the recombined data live
    opts->metafits             = NULL;  // filename of the metafits file for the target observation
    opts->coarse_chan_str      = NULL;  // Absolute or relative coarse channel
    opts->out_psrfits          = false; // Output PSRFITS format
    opts->out_vdif             = false; // Output VDIF format
    opts->out_filterbank       = false; // Output FILTERBANK format
    opts->out_fine             = false; // Output fine channelised data
    opts->out_coarse           = false; // Output coarse channelised data
    opts->analysis_filter      = NULL;
    opts->synth_filter         = NULL;
    opts->max_sec_per_file     = 200;   // Number of seconds per fits files
    //opts->gpu_mem_GB           = 0.0;
    opts->custom_flags         = NULL;
    opts->nchunks              = 1;
    opts->smart                = false;

    opts->cal_metafits         = NULL;  // filename of the metafits file for the calibration observation
    opts->caldir               = NULL;  // The path to where the calibration solutions live
    opts->cal_type             = CAL_RTS;
    opts->ref_ant              = NULL;
    opts->keep_cross_terms     = false;
    opts->phase_offset         = 0.0;
    opts->phase_slope          = 0.0;
    opts->custom_pq_correction = false;
    opts->use_bandpass         = false; // use the Bandpass calibration solutions

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"bandpass",        no_argument,       0, 'B'},
                {"out-format",      required_argument, 0, 'o'},
                {"max_t",           required_argument, 0, 't'},
                {"analysis_filter", required_argument, 0, 'A'},
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
                //{"gpu-mem",         required_argument, 0, 'g'},
                {"nchunks",         required_argument, 0, 'n'},
                {"smart",           no_argument,       0, 's'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "A:b:Bc:C:d:e:f:F:hm:n:o:OP:R:sS:t:T:U:VX",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c)
            {
                case 'A':
                    opts->analysis_filter = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->analysis_filter, optarg );
                    break;
                case 'b':
                    opts->begin_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->begin_str, optarg );
                    break;
                case 'B':
                    opts->use_bandpass = true;
                    break;
                case 'c':
                    opts->cal_metafits = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->cal_metafits, optarg );
                    break;
                case 'C':
                    opts->caldir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->caldir, optarg );
                    break;
                case 'd':
                    opts->datadir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->datadir, optarg );
                    break;
                case 'f':
                    opts->coarse_chan_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->coarse_chan_str, optarg );
                    break;
                case 'F':
                    opts->custom_flags = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->custom_flags, optarg );
                    break;
                //case 'g':
                //    opts->gpu_mem_GB = atof(optarg);
                //    break;
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'n':
                    opts->nchunks = atoi(optarg);
                    break;
                case 'o':
                    if (strcmp( optarg, "PSRFITS" ) == 0)
                    {
#ifdef HAVE_PSRFITS
                        opts->out_psrfits = true;
                        opts->out_fine = true;
#else
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "VCSBeam compiled without the PSRFITS_UTILS library, "
                                "so PSRFITS output is not available.\n" );
                        exit(EXIT_FAILURE);
#endif
                    }
                    else if (strcmp( optarg, "VDIF" ) == 0)
                    {
#ifdef HAVE_VDIF
                        opts->out_vdif = true;
                        opts->out_coarse = true;
#else
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "VCSBeam compiled without the VDIFIO library, "
                                "so VDIF output is not available.\n" );
                        exit(EXIT_FAILURE);
#endif
                    }
                    else if (strcmp( optarg, "FILTERBANK" ) == 0)
                    {
#ifdef HAVE_FILTERBANK
                        opts->out_filterbank = true;
                        opts->out_fine = true;
#else
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "VCSBeam compiled without the SIGPROC library, "
                                "so FILTERBANK output is not available.\n" );
                        exit(EXIT_FAILURE);
#endif
                    }
                    else
                    {
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "unrecognised output format \"%s\"\n", optarg );
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'O':
                    opts->cal_type = CAL_OFFRINGA;
                    break;
                case 'p':
                    opts->out_fine = true;
                    break;
                case 'P':
                    opts->pointings_file = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->pointings_file, optarg );
                    break;
                case 'R':
                    opts->ref_ant = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->ref_ant, optarg );
                    break;
                case 's':
                    opts->smart = true;
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
                    if (sscanf( optarg, "%lf,%lf", &(opts->phase_slope), &(opts->phase_offset) ) != 2)
                    {
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "cannot parse -U option (\"%s\") as \"FLOAT,FLOAT\"\n", optarg );
                        exit(EXIT_FAILURE);
                    }
                    opts->custom_pq_correction = true;
                    break;
                case 'v':
                    opts->out_coarse = true;
                    break;
                case 'V':
                    printf( "MWA Beamformer %s\n", VCSBEAM_VERSION);
                    exit(0);
                    break;
                case 'X':
                    opts->keep_cross_terms = true;
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
    assert( opts->pointings_file  != NULL );
    assert( opts->caldir          != NULL );
    assert( opts->metafits        != NULL );
    assert( opts->cal_metafits    != NULL );
    assert( opts->out_psrfits | opts->out_vdif | opts->out_filterbank );
    assert( opts->nchunks         >= 1 );

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

    if (opts->analysis_filter == NULL)
    {
        opts->analysis_filter = (char *)malloc( 8 );
        strcpy( opts->analysis_filter, "FINEPFB" );
    }

    if (opts->synth_filter == NULL)
    {
        opts->synth_filter = (char *)malloc( 6 );
        strcpy( opts->synth_filter, "LSQ12" );
    }

}



void write_step( vcsbeam_context *vm,
#ifdef HAVE_PSRFITS
        mpi_psrfits *mpfs,
#endif
#ifdef HAVE_VDIF
        struct vdifinfo *vf, vdif_header *vhdr, float *data_buffer_vdif,
#endif
        struct make_tied_array_beam_opts *opts )
{
    int p;
    for (p = 0; p < vm->npointing; p++)
    {
        logger_start_stopwatch( vm->log, "write", true );

#ifdef HAVE_PSRFITS
        if (opts->out_psrfits)
        {

            // Write out the spliced channels
            wait_splice_psrfits( &(mpfs[p]) );

            if (vm->coarse_chan_idx == mpfs[p].writer_id)
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
#endif

#ifdef HAVE_VDIF
        if (opts->out_vdif)
        {
            vdif_write_second( &vf[p], vhdr,
                    data_buffer_vdif + p * vf->sizeof_buffer );
        }
#endif

        logger_stop_stopwatch( vm->log, "write" );
    }

}
