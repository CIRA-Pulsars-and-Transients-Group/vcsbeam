/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

// Standard library
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>

// Non-standard dependencies
#include <mwalib.h>
#include <cuda_runtime.h>

// Local includes
#include "jones.h"
#include "beam_psrfits.h"
#include "metadata.h"
#include "form_beam.h"
#include "geometry.h"
#include "performance.h"

#define MAX_COMMAND_LENGTH 1024

/***********************************************
* Struct for handling the command line options *
***********************************************/

struct cmd_line_opts {
    char              *begin_str;     // Absolute or relative GPS time -- when to start beamforming
    unsigned long int  nseconds;      // How many seconds to process
    char              *datadir;       // The path to where the recombined data live
    char              *metafits;      // Filename of the metafits file
    char              *coarse_chan_str;   // Absolute or relative coarse channel number
    int                ncoarse_chans; // How many coarse channels to process
    char              *outfile;       // Base name of the output PSRFITS file
    int                max_sec_per_file;    // Number of seconds per fits file
};

/*************************************
 * Function prototypes and constants *
 *************************************/

void make_incoh_beam_parse_cmdline(
        int argc, char **argv, struct cmd_line_opts *opts );

const uintptr_t outpol_incoh = 1;  // ("I")

/********
 * MAIN *
 *******/

int main(int argc, char **argv)
{
    // Parse command line arguments
    struct cmd_line_opts opts;
    make_incoh_beam_parse_cmdline( argc, argv, &opts );

    // Start a logger for output messages and time-keeping
    logger *log = create_logger( stdout );
    logger_add_stopwatch( log, "read" );
    logger_add_stopwatch( log, "calc" );
    logger_add_stopwatch( log, "write" );
    char log_message[MAX_COMMAND_LENGTH];

    // Create an mwalib metafits context and associated metadata
    logger_timed_message( log, "Creating metafits and voltage contexts via MWALIB" );

    char error_message[ERROR_MESSAGE_LEN];

    MetafitsContext  *obs_context  = NULL;
    MetafitsMetadata *obs_metadata = NULL;
    get_mwalib_metafits_metadata( opts.metafits, &obs_metadata, &obs_context );

    unsigned long int begin_gps = parse_begin_string( obs_metadata, opts.begin_str );
    uintptr_t begin_coarse_chan_idx = parse_coarse_chan_string( obs_metadata, opts.coarse_chan_str );

    VoltageContext   *vcs_context  = NULL;
    VoltageMetadata  *vcs_metadata = NULL;
    get_mwalib_voltage_metadata( &vcs_metadata, &vcs_context, &obs_metadata, obs_context,
            begin_gps, opts.nseconds, opts.datadir, begin_coarse_chan_idx, &opts.ncoarse_chans );

    // Create some "shorthand" variables for code brevity
    uintptr_t    ntimesteps      = vcs_metadata->num_provided_timesteps;
    uintptr_t    nchans          = obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t    ninputs         = obs_metadata->num_rf_inputs;
    unsigned int nsamples        = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;
    uintptr_t    data_size       = vcs_metadata->num_voltage_blocks_per_timestep * vcs_metadata->voltage_block_size_bytes;
    uintptr_t    incoh_size      = nchans * nsamples * sizeof(float);
    uintptr_t    timestep_idx;
    uintptr_t    timestep;
    uint64_t     gps_second;
    double       ra_hours        = obs_metadata->ra_tile_pointing_deg / 15.0;
    double       dec_degs        = obs_metadata->dec_tile_pointing_deg;

    // Allocate memory
    logger_timed_message( log, "Allocate host and device memory" );

    uint8_t *data, *d_data;
    float *incoh, *d_incoh;

    allocate_input_output_arrays( (void **)&data, (void **)&d_data, data_size );
    allocate_input_output_arrays( (void **)&incoh, (void **)&d_incoh, incoh_size );

    // Get pointing geometry information
    struct beam_geom beam_geom_vals;

    double mjd = obs_metadata->sched_start_mjd;
    calc_beam_geom( ra_hours, dec_degs, mjd, &beam_geom_vals );

    struct psrfits pf;

    // COARSE CHANNEL DEPENDENT CODE BEGINS HERE
    for (uintptr_t coarse_chan_idx = begin_coarse_chan_idx; coarse_chan_idx < begin_coarse_chan_idx + opts.ncoarse_chans; coarse_chan_idx++)
    {
        // Populate the PSRFITS header struct
        sprintf( log_message, "Preparing header for output PSRFITS (receiver channel %lu)",
                obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
        logger_message( log, log_message );

        populate_psrfits_header( &pf, obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
                outpol_incoh, &beam_geom_vals, opts.outfile, false );

        // Begin the main loop: go through data one second at a time

        logger_message( log, "\n*****BEGIN BEAMFORMING*****" );

        for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
        {
            logger_message( log, "" ); // Print a blank line

            // Get the next gps second 
            timestep = vcs_metadata->provided_timestep_indices[timestep_idx];
            gps_second = vcs_metadata->timesteps[timestep].gps_time_ms / 1000;

            // Read in data from next file
            sprintf( log_message, "[%lu/%lu] Reading in data for gps second %ld",
                    timestep_idx+1, ntimesteps, gps_second );
            logger_timed_message( log, log_message );

            logger_start_stopwatch( log, "read" );
            if (mwalib_voltage_context_read_second(
                        vcs_context,
                        gps_second,
                        1,
                        coarse_chan_idx,
                        data,
                        data_size,
                        error_message,
                        ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
            {
                fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s\n", error_message );
                exit(EXIT_FAILURE);
            }
            logger_stop_stopwatch( log, "read" );

            // Form the incoherent beam
            sprintf( log_message, "[%lu/%lu] Calculating beam", timestep_idx+1, ntimesteps );
            logger_timed_message( log, log_message );

            logger_start_stopwatch( log, "calc" );

            cu_form_incoh_beam(
                    data, d_data, data_size,
                    incoh, d_incoh, incoh_size,
                    nsamples, nchans, ninputs );

            logger_stop_stopwatch( log, "calc" );


            // Write out to file
            sprintf( log_message, "[%lu/%lu] Writing data to file", timestep_idx+1, ntimesteps );
            logger_timed_message( log, log_message );

            logger_start_stopwatch( log, "write" );

            psrfits_write_second( &pf, incoh, nchans, outpol_incoh, 0 );

            logger_stop_stopwatch( log, "write" );
        }

        logger_message( log, "\n*****END BEAMFORMING*****\n" );

        // Clean up psrfits struct, ready for next round
        free_psrfits( &pf );
    }

    logger_report_all_stats( log );
    logger_message( log, "" );

    logger_timed_message( log, "Starting memory clean-up" );

    // Free up memory
    free( opts.begin_str );
    free( opts.coarse_chan_str );
    free( opts.datadir   );
    free( opts.metafits  );

    free_input_output_arrays( data, d_data );
    free_input_output_arrays( incoh, d_incoh );

    // Clean up memory associated with mwalib
    mwalib_metafits_metadata_free( obs_metadata );
    mwalib_voltage_metadata_free( vcs_metadata );
    mwalib_voltage_context_free( vcs_context );

    // Free the logger
    logger_timed_message( log, "Exiting successfully" );
    destroy_logger( log );

    return EXIT_SUCCESS;
}


void usage()
{
    printf( "\nusage: make_incoh_beam [OPTIONS]\n");

    printf( "\nREQUIRED OPTIONS\n\n"
            "\t-m, --metafits=FILE       FILE is the metafits file for the target observation\n"
           );

    printf( "\nOUTPUT OPTIONS\n\n"
            "\t-b, --begin=GPSTIME       Begin time of observation, in GPS seconds\n"
            "\t                          If GPSTIME starts with a '+' or a '-', then the time\n"
            "\t                          is taken relative to the start or end of the observation\n"
            "\t                          respectively. [default: \"+0\"]\n"
            "\t-d, --data-location=PATH  PATH is the directory containing the recombined data\n"
            "\t                          [default: current directory]\n"
            "\t-f, --coarse-chan=CHAN     Coarse channel number\n"
            "\t                           If CHAN starts with a '+' or a '-', then the channel is taken\n"
            "\t                           relative to the first or last channel in the observation\n"
            "\t                           respectively. Otherwise, it is treated as a receiver channel number\n"
            "\t                           (0-255) [default: \"+0\"]\n"
            "\t-F, --nchans=VAL           Process VAL coarse channels\n"
            "\t-o, --outfile             The base name for the output PSRFITS file\n"
            "\t                          [default: \"<PROJECT>_<OBSID>_incoh_ch<CHAN>\"]\n"
            "\t-S, --max_output_t=SECS   Maximum number of SECS per output FITS file [default: 200]\n"
            "\t-T, --nseconds=VAL        Process VAL seconds of data [default: as many as possible]\n"
           );

    printf( "\nOTHER OPTIONS\n\n"
            "\t-h, --help                Print this help and exit\n"
            "\t-V, --version             Print version number and exit\n\n"
           );
}



void make_incoh_beam_parse_cmdline(
        int argc, char **argv, struct cmd_line_opts *opts )
{
    // Set defaults for command line options
    opts->begin_str        = NULL; // Absolute or relative GPS time -- when to start beamforming
    opts->nseconds         = -1;   // How many seconds to process (-1 = as many as possible)
    opts->datadir          = NULL; // The path to where the recombined data live
    opts->metafits         = NULL; // Filename of the metafits file for the target observation
    opts->outfile          = NULL; // Base name of the output PSRFITS file
    opts->coarse_chan_str  = NULL; // Absolute or relative coarse channel
    opts->ncoarse_chans    = -1;   // How many coarse channels to process
    opts->max_sec_per_file = 200;  // Number of seconds per fits files

    if (argc > 1)
    {
        int c;
        while (1)
        {
            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"data-location",   required_argument, 0, 'd'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"nchans",          required_argument, 0, 'F'},
                {"help",            required_argument, 0, 'h'},
                {"metafits",        required_argument, 0, 'm'},
                {"outfile",         required_argument, 0, 'o'},
                {"max_output_t",    required_argument, 0, 'S'},
                {"nseconds",        required_argument, 0, 'T'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "b:d:f:F:hm:o:S:T:V",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'b':
                    opts->begin_str = strdup(optarg);
                    break;
                case 'd':
                    opts->datadir = strdup(optarg);
                    break;
                case 'f':
                    opts->coarse_chan_str = strdup(optarg);
                    break;
                case 'F':
                    opts->ncoarse_chans = atoi(optarg);
                    break;
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'o':
                    opts->outfile = strdup(optarg);
                    break;
                case 'S':
                    opts->max_sec_per_file = atoi(optarg);
                    break;
                case 'T':
                    opts->nseconds = atol(optarg);
                    if (opts->nseconds <= 0)
                    {
                        fprintf( stderr, "error: make_incoh_beam_parse_cmdline: "
                                "-%c argument must be >= 1\n", c );
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'V':
                    printf( "MWA Beamformer %s\n", VERSION_BEAMFORMER );
                    exit(EXIT_SUCCESS);
                    break;
                default:
                    fprintf( stderr, "error: make_incoh_beam_parse_cmdline: "
                                    "unrecognised option '%s'\n", optarg );
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
    assert( opts->metafits     != NULL );

    if (opts->datadir == NULL)
        opts->datadir = strdup( "." );

    if (opts->begin_str == NULL)
        opts->begin_str = strdup( "+0" );

    if (opts->coarse_chan_str == NULL)
        opts->coarse_chan_str = strdup( "+0" );
}


