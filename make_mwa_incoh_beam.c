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
#include "metadata.h"
#include "psrfits.h"
#include "form_beam.h"
#include "geometric_delay.h"

#include <cuda_runtime.h>
#include "ipfb.h"
#include "performance.h"

#define MAX_COMMAND_LENGTH 1024

/***********************************************
* Struct for handling the command line options *
***********************************************/

struct cmd_line_opts {
    // Variables for required options
    unsigned long int  begin;         // GPS time -- when to start beamforming
    unsigned long int  end;           // GPS time -- when to stop beamforming
    char              *datadir;       // The path to where the recombined data live
    char              *metafits;      // filename of the metafits file
    uintptr_t          rec_channel;   // 0 - 255 receiver 1.28MHz channel

    // Other options
    int                max_sec_per_file;    // Number of seconds per fits files
};

/*************************************
 * Function prototypes and constants *
 *************************************/

void make_incoh_beam_parse_cmdline(
        int argc, char **argv, struct cmd_line_opts *opts );

void allocate_input_output_arrays( void **data, void **d_data, size_t size );
void free_input_output_arrays( void *data, void *d_data );

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
    MetafitsMetadata *obs_metadata = NULL;
    VoltageMetadata  *vcs_metadata = NULL;
    VoltageContext   *vcs_context  = NULL;

    get_mwalib_metadata( &obs_metadata, &vcs_metadata, &vcs_context, NULL,
            opts.metafits, NULL, opts.begin, opts.end, opts.datadir, opts.rec_channel );

    // Create some "shorthand" variables for code brevity
    uintptr_t    ntimesteps      = vcs_metadata->num_common_timesteps;
    uintptr_t    nchans          = obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t    ninputs         = obs_metadata->num_rf_inputs;
    unsigned int nsamples        = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;
    int          coarse_chan_idx = 0; /* Value is fixed for now (i.e. each call of make_beam only
                                     ever processes one coarse chan. However, in the future,
                                     this should be flexible, with mpi or threads managing
                                     different coarse channels. */
    int          coarse_chan     = vcs_metadata->common_coarse_chan_indices[coarse_chan_idx];
    uintptr_t    data_size       = vcs_metadata->num_voltage_blocks_per_timestep * vcs_metadata->voltage_block_size_bytes;
    uintptr_t    incoh_size      = nchans * nsamples * sizeof(float);
    uintptr_t    timestep_idx;
    uintptr_t    timestep;
    uint64_t     gps_second;
    double       ra_hours        = obs_metadata->ra_tile_pointing_deg / 15.0;
    double       dec_degs        = obs_metadata->dec_tile_pointing_deg;

    // Populate the PSRFITS header struct
    logger_timed_message( log, "Preparing header for output PSRFITS" );

    struct beam_geom beam_geom_vals;

    double mjd = obs_metadata->sched_start_mjd;
    calc_beam_geom( ra_hours, dec_degs, mjd, &beam_geom_vals );

    struct psrfits pf;

    populate_psrfits_header( &pf, obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
            outpol_incoh, &beam_geom_vals, false );

    // Allocate memory
    logger_timed_message( log, "Allocate host and device memory" );

    // Raw data:
    uint8_t *data, *d_data;
    float *incoh, *d_incoh;

    allocate_input_output_arrays( (void **)&data, (void **)&d_data, data_size );
    allocate_input_output_arrays( (void **)&incoh, (void **)&d_incoh, incoh_size );

    // Begin the main loop: go through data one second at a time

    logger_message( log, "\n*****BEGIN BEAMFORMING*****" );

    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        logger_message( log, "" ); // Print a blank line

        // Get the next gps second 
        timestep = vcs_metadata->common_timestep_indices[timestep_idx];
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
                    coarse_chan,
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

        psrfits_write_second( &pf, incoh,
                nchans, outpol_incoh, 0 );

        logger_stop_stopwatch( log, "write" );
    }

    logger_message( log, "\n*****END BEAMFORMING*****\n" );
    logger_report_all_stats( log );
    logger_message( log, "" );

    logger_timed_message( log, "Starting memory clean-up" );

    // Free up memory
    free_input_output_arrays( data, d_data );
    free_input_output_arrays( incoh, d_incoh );

    free( opts.datadir        );
    free( opts.metafits       );

    free_psrfits( &pf );

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
            "\t-b, --begin=GPSTIME       Begin time of observation, in GPS seconds\n"
            "\t-e, --end=GPSTIME         End time of observation, in GPS seconds\n"
            "\t-d, --data-location=PATH  PATH is the directory containing the recombined data\n"
            "\t-m, --metafits=FILE       FILE is the metafits file for the target observation\n"
            "\t-f, --coarse-chan=N       Receiver coarse channel number (0-255)\n"
           );

    printf( "\nOUTPUT OPTIONS\n\n"
            "\t-t, --max_t               Maximum number of seconds per output FITS file [default: 200]\n"
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
    opts->begin            = 0;    // GPS time -- when to start beamforming
    opts->end              = 0;    // GPS time -- when to stop beamforming
    opts->datadir          = NULL; // The path to where the recombined data live
    opts->metafits         = NULL; // filename of the metafits file for the target observation
    opts->rec_channel      = -1;   // 0 - 255 receiver 1.28MHz channel
    opts->max_sec_per_file = 200;  // Number of seconds per fits files

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"data-location",   required_argument, 0, 'd'},
                {"end",             required_argument, 0, 'e'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"help",            required_argument, 0, 'h'},
                {"metafits",        required_argument, 0, 'm'},
                {"max_t",           required_argument, 0, 't'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "b:d:e:f:hm:t:V",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'b':
                    opts->begin = atol(optarg);
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
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 't':
                    opts->max_sec_per_file = atoi(optarg);
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
    assert( opts->begin        != 0    );
    assert( opts->end          != 0    );
    assert( opts->datadir      != NULL );
    assert( opts->metafits     != NULL );
}


void allocate_input_output_arrays( void **data, void **d_data, size_t size )
{
    cudaMallocHost( data, size );
    cudaCheckErrors( "cudaMallocHost() failed" );

    cudaMalloc( d_data, size );
    cudaCheckErrors( "cudaMalloc() failed" );
}

void free_input_output_arrays( void *data, void *d_data )
{
    cudaFreeHost( data );
    cudaFree( d_data );
}
