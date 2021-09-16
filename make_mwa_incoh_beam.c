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
#include "vcsbeam.h"
#include "jones.h"
#include "beam_psrfits.h"
#include "metadata.h"
#include "form_beam.h"
#include "geometry.h"

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

void read_step( VoltageContext *vcs_context, uint64_t gps_second, uintptr_t
        coarse_chan_idx, uint8_t *data, uintptr_t data_size, logger *log );
void write_step( mpi_psrfits *mpf, logger *log );

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
    const int chans_per_proc = 1;

    // Parse command line arguments
    struct cmd_line_opts opts;
    make_incoh_beam_parse_cmdline( argc, argv, &opts );

    // Start a logger for output messages and time-keeping
    logger *log = create_logger( stdout, mpi_proc_id );
    logger_add_stopwatch( log, "read", "Reading in data" );
    logger_add_stopwatch( log, "calc", "Calculating the incoherent beam" );
    logger_add_stopwatch( log, "splice", "Splicing coarse channels together" );
    if( mpi_proc_id == writer )
        logger_add_stopwatch( log, "write", "Writing data to file" );
    char log_message[MAX_COMMAND_LENGTH];

    // Create an mwalib metafits context and associated metadata
    logger_timed_message( log, "Creating metafits and voltage contexts via MWALIB" );

    vcsbeam_metadata *vm = init_vcsbeam_metadata(
        opts.metafits, NULL,
        opts.coarse_chan_str, chans_per_proc, mpi_proc_id,
        opts.begin_str, opts.nseconds, 0,
        opts.datadir );

    // Create some "shorthand" variables for code brevity
    uintptr_t    nchans          = vm->obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t    ninputs         = vm->obs_metadata->num_rf_inputs;
    unsigned int nsamples        = vm->sample_rate;
    uintptr_t    data_size       = vm->bytes_per_second;
    uintptr_t    incoh_size      = nchans * nsamples * sizeof(float);
    uintptr_t    Iscaled_size    = nchans * nsamples * sizeof(uint8_t);
    double       ra_hours        = vm->obs_metadata->ra_tile_pointing_deg / 15.0;
    double       dec_degs        = vm->obs_metadata->dec_tile_pointing_deg;
    uintptr_t    timestep_idx;
    uint64_t     gps_second;

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

    double mjd = vm->obs_metadata->sched_start_mjd;
    calc_beam_geom( ra_hours, dec_degs, mjd, &beam_geom_vals );

    mpi_psrfits mpf;
    init_mpi_psrfits(
        &mpf,
        vm->obs_metadata,
        vm->vcs_metadata,
        world_size,
        mpi_proc_id,
        opts.max_sec_per_file,
        outpol_incoh,
        &beam_geom_vals,
        opts.outfile,
        writer,
        false );

    // Begin the main loop: go through data one second at a time

    logger_message( log, "\n*****BEGIN BEAMFORMING*****" );

    uintptr_t ntimesteps = vm->num_gps_seconds_to_process;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        // Get the next gps second 
        gps_second = vm->gps_seconds_to_process[timestep_idx];

        sprintf( log_message, "---Processing GPS second %ld [%lu/%lu]---",
                gps_second, timestep_idx+1, ntimesteps );
        logger_message( log, log_message );

        // Read in data from next file
        read_step( vm->vcs_context, gps_second, vm->coarse_chan_idxs_to_process[0], data, data_size, log );

        // The writing (of the previous second) is put here in order to
        // allow the possibility that it can overlap with the reading step.
        // Because of this, another "write" has to happen after this loop
        // has terminated
        if (timestep_idx > 0) // i.e. don't do this the first time around
        {
            write_step( &mpf, log );
        }

        // Form the incoherent beam
        logger_start_stopwatch( log, "calc", true );

        cu_form_incoh_beam(
                data, d_data, data_size,
                d_incoh,
                nsamples, nchans, ninputs,
                mpf.coarse_chan_pf.sub.dat_offsets, d_offsets,
                mpf.coarse_chan_pf.sub.dat_scales, d_scales,
                mpf.coarse_chan_pf.sub.data, d_Iscaled, Iscaled_size
                );

        logger_stop_stopwatch( log, "calc" );


        // Splice the channels together
        logger_start_stopwatch( log, "splice", true );

        gather_splice_psrfits( &mpf );

        logger_stop_stopwatch( log, "splice" );
    }

    // Write out the last second's worth of data
    write_step( &mpf, log );

    logger_message( log, "\n*****END BEAMFORMING*****\n" );

    logger_report_all_stats( log );
    logger_message( log, "" );

    logger_timed_message( log, "Starting memory clean-up" );

    // Free up memory
    free_mpi_psrfits( &mpf );

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
    destroy_vcsbeam_metadata( vm );

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

void read_step( VoltageContext *vcs_context, uint64_t gps_second,
        uintptr_t coarse_chan_idx, uint8_t *data, uintptr_t data_size,
        logger *log )
{
    // Read second's worth of data from file
    logger_start_stopwatch( log, "read", true );

    char error_message[ERROR_MESSAGE_LEN];
    if (mwalib_voltage_context_read_second(
                vcs_context,
                gps_second,
                1,
                coarse_chan_idx,
                data,
                data_size,
                error_message,
                ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    logger_stop_stopwatch( log, "read" );
}


void write_step( mpi_psrfits *mpf, logger *log )
{
    wait_splice_psrfits( mpf );

    int mpi_proc_id;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_proc_id );

    if (mpi_proc_id == mpf->writer_id)
    {
        // Write second's worth of data to file
        logger_start_stopwatch( log, "write", true );

        if (psrfits_write_subint( &(mpf->spliced_pf) ) != 0)
        {
            fprintf(stderr, "error: Write PSRFITS subint failed. File exists?\n");
            exit(EXIT_FAILURE);
        }
        mpf->spliced_pf.sub.offs = roundf(mpf->spliced_pf.tot_rows * mpf->spliced_pf.sub.tsubint) + 0.5*mpf->spliced_pf.sub.tsubint;
        mpf->spliced_pf.sub.lst += mpf->spliced_pf.sub.tsubint;

        logger_stop_stopwatch( log, "write" );
    }
}
