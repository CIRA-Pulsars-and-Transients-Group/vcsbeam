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

// Local includes
#include "vcsbeam.h"

struct fine_pfb_offline_opts {
    char              *begin_str;        // Absolute or relative GPS time -- when to start beamforming
    unsigned long int  nseconds;         // How many seconds to process
    char              *datadir;          // The path to where the recombined data live
    char              *metafits;         // filename of the metafits file
    char              *coarse_chan_str;  // Absolute or relative coarse channel number
    char              *analysis_filter;  // Which analysis filter to use
    int                nchunks;          // Split each second into this many processing chunks
};

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/

void usage();
void fine_pfb_offline_parse_cmdline( int argc, char **argv, struct fine_pfb_offline_opts *opts );

/********
 * MAIN *
 ********/

int main( int argc, char *argv[] )
{
    // Parse command line arguments
    struct fine_pfb_offline_opts opts;
    fine_pfb_offline_parse_cmdline( argc, argv, &opts );

    // Set up the VCS metadata struct
    bool use_mpi = false;
    vcsbeam_context *vm = vmInit( use_mpi );

    vmLoadObsMetafits( vm, opts.metafits );
    vmBindObsData( vm,
        opts.coarse_chan_str, 1, 0,
        opts.begin_str, opts.nseconds, 0,
        opts.datadir );

    vmPrintTitle( vm, "Offline PFB" );

    // This only works for new-style (coarse-channelised) MWAX data
    if (vm->obs_metadata->mwa_version != VCSMWAXv2)
    {
        fprintf( stderr, "error: observation %u does not appear to be "
                "coarse-channelised MWAX data\n",
                vm->obs_metadata->obs_id );
        exit(EXIT_FAILURE);
    }

    // Load the filter
    int K = 128; // The number of desired output channels
    vmLoadFilter( vm, opts.analysis_filter, ANALYSIS_FILTER, K );

    // Create and init the PFB struct
    int M = K; // The filter stride (M = K <=> "critically sampled PFB")
    vm->chunks_per_second = opts.nchunks;
    vmInitForwardPFB( vm, M, PFB_SMART | PFB_MALLOC_HOST_OUTPUT );

    // Let's try it out on one second of data
    char filename[128];
    vm_error status;
    while ((status = vmReadNextSecond( vm )) == VM_SUCCESS)
    {
        // Actually do the PFB
        vmExecuteForwardPFB( vm );

        // Write out the answer to a file
        vmWritePFBOutputToFile( vm );
    }

    // Report performance statistics
    vmReportPerformanceStats( vm );

    // Free memory
    destroy_vcsbeam_context( vm );

    // Exit gracefully
    return EXIT_SUCCESS;
}

/************************
 * FUNCTION DEFINITIONS *
 ************************/

void usage()
{
    printf( "\nusage: fine_pfb_offline [OPTIONS]\n");

    printf( "\nREQUIRED OPTIONS\n\n"
            "\t-m, --metafits=FILE        FILE is the metafits file for the target observation\n"
          );

    printf( "\nOPTIONAL OPTIONS\n\n"
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
            "\t-A, --analysis_filter=FILTER  Apply the named filter during fine channelisation.\n"
            "\t                           File [RUNTIME_DIR]/FILTER.dat must exist [default: FINEPFB]\n"
            "\t-T, --nseconds=VAL         Process VAL seconds of data [default: as many as possible]\n"
            "\t-n, --nchunks=VAL          Split each second's worth of data into VAL processing chunks\n"
            "\t                           [default: 1]\n"
          );

    printf( "\nOTHER OPTIONS\n\n"
            "\t-h, --help                 Print this help and exit\n"
            "\t-V, --version              Print version number and exit\n\n"
          );
}

void fine_pfb_offline_parse_cmdline( int argc, char **argv, struct fine_pfb_offline_opts *opts )
{
    // Set defaults
    opts->begin_str          = NULL;  // Absolute or relative GPS time -- when to start beamforming
    opts->nseconds           = -1;    // How many seconds to process (-1 = as many as possible)
    opts->datadir            = NULL;  // The path to where the recombined data live
    opts->metafits           = NULL;  // filename of the metafits file for the target observation
    opts->coarse_chan_str    = NULL;  // Absolute or relative coarse channel
    opts->analysis_filter    = NULL;
    opts->nchunks            = 1;

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"analysis_filter", required_argument, 0, 'A'},
                {"begin",           required_argument, 0, 'b'},
                {"data-location",   required_argument, 0, 'd'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"help",            required_argument, 0, 'h'},
                {"metafits",        required_argument, 0, 'm'},
                {"nchunks",         required_argument, 0, 'n'},
                {"nseconds",        required_argument, 0, 'T'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "A:b:d:f:hm:n:T:V",
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
                case 'd':
                    opts->datadir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->datadir, optarg );
                    break;
                case 'f':
                    opts->coarse_chan_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->coarse_chan_str, optarg );
                    break;
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'm':
                    opts->metafits = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->metafits, optarg );
                    break;
                case 'n':
                    opts->nchunks = atoi(optarg);
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
                case 'V':
                    printf( "MWA Beamformer %s\n", VCSBEAM_VERSION);
                    exit(0);
                    break;
                default:
                    fprintf( stderr, "error: fine_pfb_offline_parse_cmdline: "
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
    // and/or set defaults
    assert( opts->metafits != NULL );
    assert( opts->nchunks >= 1 );

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
}
