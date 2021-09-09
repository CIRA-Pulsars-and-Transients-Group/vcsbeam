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

// Local includes
#include "pfb.h"
#include "metadata.h"
#include "performance.h"

struct fine_pfb_offline_opts {
    char              *begin_str;        // Absolute or relative GPS time -- when to start beamforming
    unsigned long int  nseconds;         // How many seconds to process
    char              *datadir;          // The path to where the recombined data live
    char              *metafits;         // filename of the metafits file
    char              *coarse_chan_str;  // Absolute or relative coarse channel number
    char              *synth_filter;     // Which synthesis filter to use
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

    // Start a logger for output messages and time-keeping
    logger *log = create_logger( stdout, PERFORMANCE_NO_MPI );
    logger_add_stopwatch( log, "read", "Reading in data" );
    logger_add_stopwatch( log, "pfb", "Performing the PFB" );
    logger_add_stopwatch( log, "write", "Writing out data to file" );
    //char log_message[MAX_COMMAND_LENGTH];

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
            "\t-S, --synth_filter=filter  Apply the named filter during high-time resolution synthesis.\n"
            "\t                           filter can be MIRROR or LSQ12.\n"
            "\t                           [default: LSQ12]\n"
            "\t-T, --nseconds=VAL         Process VAL seconds of data [default: as many as possible]\n"
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
    opts->synth_filter       = NULL;

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"data-location",   required_argument, 0, 'd'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"help",            required_argument, 0, 'h'},
                {"metafits",        required_argument, 0, 'm'},
                {"synth_filter",    required_argument, 0, 'S'},
                {"nseconds",        required_argument, 0, 'T'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "b:d:f:hm:S:T:V",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c)
            {
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
                    opts->metafits = strdup(optarg);
                    break;
                case 'S':
                    opts->synth_filter = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->synth_filter, optarg );
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
                    printf( "MWA Beamformer %s\n", VERSION_BEAMFORMER);
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
