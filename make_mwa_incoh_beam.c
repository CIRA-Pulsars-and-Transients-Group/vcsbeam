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
#include <mpi.h>

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
    // Initialise MPI
    MPI_Init( NULL, NULL );
    int world_size, mpi_proc_id;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_id);

    const int writer = 0; // Designate process 0 to write out the files
    int ncoarse_chans = world_size; // Each process handles one coarse channel

    // Parse command line arguments
    struct cmd_line_opts opts;
    make_incoh_beam_parse_cmdline( argc, argv, &opts );

    // Start a logger for output messages and time-keeping
    logger *log = create_logger( stdout, mpi_proc_id );
    logger_add_stopwatch( log, "read" );
    logger_add_stopwatch( log, "calc" );
    logger_add_stopwatch( log, "splice" );
    if( mpi_proc_id == writer )
        logger_add_stopwatch( log, "write" );
    char log_message[MAX_COMMAND_LENGTH];

    // Create an mwalib metafits context and associated metadata
    logger_timed_message( log, "Creating metafits and voltage contexts via MWALIB" );

    char error_message[ERROR_MESSAGE_LEN];

    MetafitsContext  *obs_context  = NULL;
    MetafitsMetadata *obs_metadata = NULL;
    get_mwalib_metafits_metadata( opts.metafits, &obs_metadata, &obs_context );

    unsigned long int begin_gps = parse_begin_string( obs_metadata, opts.begin_str );

    uintptr_t coarse_chan_idx = parse_coarse_chan_string( obs_metadata, opts.coarse_chan_str ) + mpi_proc_id;

    sprintf( log_message, "rank = %d, coarse_chan_idx = %lu\n", mpi_proc_id, coarse_chan_idx );
    logger_timed_message( log, log_message );

    VoltageContext   *vcs_context  = NULL;
    VoltageMetadata  *vcs_metadata = NULL;
    get_mwalib_voltage_metadata( &vcs_metadata, &vcs_context, &obs_metadata, obs_context,
            begin_gps, opts.nseconds, opts.datadir, coarse_chan_idx, 1 );

    // Create some "shorthand" variables for code brevity
    uintptr_t    ntimesteps      = vcs_metadata->num_provided_timesteps;
    uintptr_t    nchans          = obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t    ninputs         = obs_metadata->num_rf_inputs;
    unsigned int nsamples        = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;
    uintptr_t    data_size       = vcs_metadata->num_voltage_blocks_per_timestep * vcs_metadata->voltage_block_size_bytes;
    uintptr_t    incoh_size      = nchans * nsamples * sizeof(float);
    uintptr_t    Iscaled_size    = nchans * nsamples * sizeof(uint8_t);
    uintptr_t    timestep_idx;
    uintptr_t    timestep;
    uint64_t     gps_second;
    double       ra_hours        = obs_metadata->ra_tile_pointing_deg / 15.0;
    double       dec_degs        = obs_metadata->dec_tile_pointing_deg;

    // Allocate memory
    logger_timed_message( log, "Allocate host and device memory" );

    uint8_t *data, *d_data;
    float *d_incoh;
    float *d_offsets;
    float *d_scales;
    uint8_t *d_Iscaled;

    allocate_input_output_arrays( (void **)&data, (void **)&d_data, data_size );

    cudaMalloc( (void **)&d_incoh, incoh_size );
    cudaCheckErrors( "cudaMalloc(d_incoh) failed" );
    cudaMalloc( (void **)&d_offsets, nchans*sizeof(float) );
    cudaCheckErrors( "cudaMalloc(d_offsets) failed" );
    cudaMalloc( (void **)&d_scales,  nchans*sizeof(float) );
    cudaCheckErrors( "cudaMalloc(d_scales) failed" );
    cudaMalloc( (void **)&d_Iscaled, Iscaled_size );
    cudaCheckErrors( "cudaMalloc(Iscaled) failed" );

    // Get pointing geometry information
    struct beam_geom beam_geom_vals;

    double mjd = obs_metadata->sched_start_mjd;
    calc_beam_geom( ra_hours, dec_degs, mjd, &beam_geom_vals );

    struct psrfits pf;
    uint8_t *spliced_buffer = NULL;

    struct psrfits spliced_pf;
    if (mpi_proc_id == writer)
    {
        logger_timed_message( log, "Preparing header for spliced PSRFITS" );
        int coarse_chan_idxs[ncoarse_chans];

        int i;
        for (i = 0; i < ncoarse_chans; i++)
            coarse_chan_idxs[i] = parse_coarse_chan_string( obs_metadata, opts.coarse_chan_str ) + i;

        populate_spliced_psrfits_header( &spliced_pf, obs_metadata, vcs_metadata, coarse_chan_idxs, ncoarse_chans,
                opts.max_sec_per_file, outpol_incoh, &beam_geom_vals, opts.outfile, false );

        spliced_buffer = (uint8_t *)malloc( spliced_pf.sub.bytes_per_subint );
    }

    // Populate the PSRFITS header struct
    sprintf( log_message, "Preparing header for output PSRFITS (receiver channel %lu)",
            obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
    logger_timed_message( log, log_message );

    populate_psrfits_header( &pf, obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
            outpol_incoh, &beam_geom_vals, opts.outfile, false );

    // Begin the main loop: go through data one second at a time

    // Create MPI vector types designed to splice the coarse channels together
    // correctly during MPI_Gather

    MPI_Datatype coarse_chan_spectrum;
    MPI_Type_contiguous( pf.hdr.nchan, MPI_BYTE, &coarse_chan_spectrum );
    MPI_Type_commit( &coarse_chan_spectrum );

    MPI_Datatype total_spectrum_type;
    MPI_Type_vector( pf.hdr.nsblk, 1, ncoarse_chans, coarse_chan_spectrum, &total_spectrum_type );
    MPI_Type_commit( &total_spectrum_type );

    MPI_Datatype spliced_type;
    MPI_Type_create_resized( total_spectrum_type, 0, pf.hdr.nchan, &spliced_type );
    MPI_Type_commit( &spliced_type );

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
                d_incoh,
                nsamples, nchans, ninputs,
                pf.sub.dat_offsets, d_offsets,
                pf.sub.dat_scales, d_scales,
                pf.sub.data, d_Iscaled, Iscaled_size
                );

        logger_stop_stopwatch( log, "calc" );


        // Write out to file

        sprintf( log_message, "[%lu/%lu] Writing data to file", timestep_idx+1, ntimesteps );
        logger_timed_message( log, log_message );

        logger_start_stopwatch( log, "splice" );

        MPI_Gather( pf.sub.data,
                pf.hdr.nsblk,
                coarse_chan_spectrum,
                spliced_pf.sub.data,
                1,
                spliced_type,
                0,
                MPI_COMM_WORLD );

        if (mpi_proc_id == writer)
        {
            // Write it to file
            logger_start_stopwatch( log, "write" );

            if (psrfits_write_subint( &spliced_pf ) != 0)
            {
                fprintf(stderr, "error: Write spliced subint failed. File exists?\n");
                exit(EXIT_FAILURE);
            }
            spliced_pf.sub.offs = roundf(spliced_pf.tot_rows * spliced_pf.sub.tsubint) + 0.5*spliced_pf.sub.tsubint;
            spliced_pf.sub.lst += spliced_pf.sub.tsubint;

            logger_stop_stopwatch( log, "write" );
        }

        logger_stop_stopwatch( log, "splice" );

    }

    logger_message( log, "\n*****END BEAMFORMING*****\n" );

    MPI_Type_free( &spliced_type );
    MPI_Type_free( &total_spectrum_type );
    MPI_Type_free( &coarse_chan_spectrum );

    // Clean up psrfits struct, ready for next round
    free_psrfits( &pf );

    logger_report_all_stats( log );
    logger_message( log, "" );

    logger_timed_message( log, "Starting memory clean-up" );

    // Free up memory
    if (mpi_proc_id == writer)
    {
        free( spliced_buffer );
    }

    free( opts.begin_str );
    free( opts.coarse_chan_str );
    free( opts.datadir   );
    free( opts.metafits  );

    free_input_output_arrays( data, d_data );

    cudaFree( d_incoh );
    cudaCheckErrors( "cudaFree(d_incoh) failed" );
    cudaFree( d_offsets );
    cudaCheckErrors( "cudaFree(d_offsets) failed" );
    cudaFree( d_scales );
    cudaCheckErrors( "cudaFree(d_scales) failed" );
    cudaFree( d_Iscaled );
    cudaCheckErrors( "cudaFree(d_Iscaled) failed" );

    // Clean up memory associated with mwalib
    mwalib_metafits_metadata_free( obs_metadata );
    mwalib_voltage_metadata_free( vcs_metadata );
    mwalib_voltage_context_free( vcs_context );

    // Free the logger
    logger_timed_message( log, "Exiting successfully" );
    destroy_logger( log );

    // Finalise MPI
    MPI_Finalize();

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
            "\t-f, --coarse-chan=CHAN    Coarse channel number\n"
            "\t                          If CHAN starts with a '+' or a '-', then the channel is taken\n"
            "\t                          relative to the first or last channel in the observation\n"
            "\t                          respectively. Otherwise, it is treated as a receiver channel number\n"
            "\t                          (0-255) [default: \"+0\"]\n"
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
                {"help",            no_argument,       0, 'h'},
                {"metafits",        required_argument, 0, 'm'},
                {"outfile",         required_argument, 0, 'o'},
                {"max_output_t",    required_argument, 0, 'S'},
                {"nseconds",        required_argument, 0, 'T'},
                {"version",         no_argument,       0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "b:d:f:hm:o:S:T:V",
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


