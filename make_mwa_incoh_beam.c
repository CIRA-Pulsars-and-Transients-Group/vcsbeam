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

double now(){
  struct timespec t;
  clock_gettime(CLOCK_REALTIME,&t);
  return (double)t.tv_sec + (double)t.tv_nsec/1000000000L;
}

#define NOW now()
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

/***********************
 * Function prototypes *
 **********************/

void make_incoh_beam_parse_cmdline(
        int argc, char **argv, struct cmd_line_opts *opts );

/********
 * MAIN *
 *******/

int main(int argc, char **argv)
{
    // A place to hold the beamformer settings
    struct cmd_line_opts opts;

    const uintptr_t outpol_incoh = 1;  // ("I")
    FILE *log = stdout; // Where the output messages go

    // Parse command line arguments
    make_incoh_beam_parse_cmdline( argc, argv, &opts );

    double begintime = NOW;

    char error_message[ERROR_MESSAGE_LEN];

    // Create an mwalib metafits context and associated metadata
    fprintf( log, "[%f]  Creating metafits and voltage contexts via MWALIB\n", NOW-begintime );

    MetafitsMetadata *obs_metadata = NULL;
    VoltageMetadata  *vcs_metadata = NULL;
    VoltageContext   *vcs_context  = NULL;

    get_mwalib_metadata( &obs_metadata, &vcs_metadata, &vcs_context, NULL,
            opts.metafits, NULL, opts.begin, opts.end, opts.datadir, opts.rec_channel );

    uintptr_t ntimesteps = vcs_metadata->num_common_timesteps;

    // Create some "shorthand" variables for code brevity
    uintptr_t nchan          = obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t ninput         = obs_metadata->num_rf_inputs;
    unsigned int sample_rate = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;


    int coarse_chan_idx = 0; /* Value is fixed for now (i.e. each call of make_beam only
                                ever processes one coarse chan. However, in the future,
                                this should be flexible, with mpi or threads managing
                                different coarse channels. */
    int coarse_chan = vcs_metadata->provided_coarse_chan_indices[coarse_chan_idx];

    fprintf( log, "[%f]  Setting up output header information\n", NOW-begintime);

    int npointing = 1;
    struct beam_geom beam_geom_vals[npointing];
    double ras_hours[npointing], decs_degs[npointing];
    ras_hours[0] = obs_metadata->ra_tile_pointing_deg / 15.0;
    decs_degs[0] = obs_metadata->dec_tile_pointing_deg;

    double mjd = obs_metadata->sched_start_mjd;
    calc_beam_geom( ras_hours, decs_degs, npointing, mjd, beam_geom_vals );

    // Populate the PSRFITS header struct
    struct psrfits  *pf_incoh;
    pf_incoh = (struct psrfits *)malloc(1 * sizeof(struct psrfits));

    populate_psrfits_header( pf_incoh, obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
            outpol_incoh, beam_geom_vals, 1, false );

    // Allocate memory for the raw data on both host and device
    long int data_size = vcs_metadata->num_voltage_blocks_per_timestep * vcs_metadata->voltage_block_size_bytes;

    uint8_t *data, *d_data;

    cudaMallocHost( (void **)&data, data_size );
    cudaCheckErrors( "cudaMallocHost(data) failed" );

    cudaMalloc( (void **)&d_data, data_size );
    cudaCheckErrors( "cudaMalloc(d_data) failed" );

    // Create output buffer arrays
    float *incoh, *d_incoh;
    size_t incoh_size = nchan * outpol_incoh * pf_incoh[0].hdr.nsblk * sizeof(float);

    cudaMallocHost( (void **)&incoh, incoh_size );
    cudaCheckErrors( "cudaMallocHost(incoh) failed" );

    cudaMalloc( (void **)&d_incoh, incoh_size );
    cudaCheckErrors( "cudaMalloc(d_incoh) failed" );

    fprintf( log, "\n[%f]  *****BEGINNING BEAMFORMING*****\n", NOW-begintime );

    // Set up timing for each section
    long read_time[ntimesteps], calc_time[ntimesteps], write_time[ntimesteps];
    uintptr_t timestep_idx;
    long timestep;
    uint64_t gps_second;

    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        timestep = vcs_metadata->common_timestep_indices[timestep_idx];
        gps_second = vcs_metadata->timesteps[timestep].gps_time_ms / 1000;

        // Read in data from next file
        clock_t start = clock();
        fprintf( log, "\n[%f] [%lu/%lu] Reading in data for gps second %ld \n", NOW-begintime,
                timestep_idx+1, ntimesteps, gps_second );

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
        read_time[timestep_idx] = clock() - start;

        fprintf( log, "[%f] [%lu/%lu] Calculating beam\n", NOW-begintime,
                                timestep_idx+1, ntimesteps);
        start = clock();

        cu_form_incoh_beam(
                data, d_data, data_size,
                incoh, d_incoh, incoh_size,
                sample_rate, nchan, ninput );

        calc_time[timestep_idx] = clock() - start;


        // Write out to file
        start = clock();
        fprintf( log, "[%f] [%lu/%lu] Writing data to file(s)\n",
                NOW-begintime, timestep_idx+1, ntimesteps );

        psrfits_write_second( &pf_incoh[0], incoh,
                nchan, outpol_incoh, 0 );

        write_time[timestep_idx] = clock() - start;
    }

    // Calculate total processing times
    float read_sum = 0, calc_sum = 0, write_sum = 0;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        read_sum  += (float) read_time[timestep_idx];
        calc_sum  += (float) calc_time[timestep_idx];
        write_sum += (float) write_time[timestep_idx];
    }
    float read_mean, calc_mean, write_mean;
    read_mean  = read_sum  / ntimesteps;
    calc_mean  = calc_sum  / ntimesteps;
    write_mean = write_sum / ntimesteps;

    // Calculate the standard deviations
    float read_std = 0, calc_std = 0, write_std = 0;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        read_std  += pow((float)read_time[timestep_idx]  - read_mean,  2);
        calc_std  += pow((float)calc_time[timestep_idx]  - calc_mean / npointing,  2);
        write_std += pow((float)write_time[timestep_idx] - write_mean / npointing, 2);
    }
    read_std  = sqrt( read_std  / ntimesteps );
    calc_std  = sqrt( calc_std  / ntimesteps / npointing );
    write_std = sqrt( write_std / ntimesteps / npointing );

    fprintf( log, "\n[%f]  *****FINISHED BEAMFORMING*****\n\n", NOW-begintime );

    fprintf( log, "[%f]  Total read  processing time: %9.3f s\n",
                     NOW-begintime, read_sum / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Mean  read  processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, read_mean / CLOCKS_PER_SEC, read_std / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Total calc  processing time: %9.3f s\n",
                     NOW-begintime, calc_sum / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Mean  calc  processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, calc_mean / CLOCKS_PER_SEC, calc_std / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Total write processing time: %9.3f s\n",
                     NOW-begintime, write_sum  * npointing / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Mean  write processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, write_mean / CLOCKS_PER_SEC, write_std / CLOCKS_PER_SEC);


    fprintf( log, "[%f]  Starting clean-up\n", NOW-begintime);

    // Free up memory
    cudaFreeHost( incoh );
    cudaFreeHost( data );

    cudaFree( d_incoh );
    cudaFree( d_data );

    free( opts.datadir        );
    free( opts.metafits       );

    free( pf_incoh[0].sub.data        );
    free( pf_incoh[0].sub.dat_freqs   );
    free( pf_incoh[0].sub.dat_weights );
    free( pf_incoh[0].sub.dat_offsets );
    free( pf_incoh[0].sub.dat_scales  );

    // Clean up memory associated with mwalib
    mwalib_metafits_metadata_free( obs_metadata );
    mwalib_voltage_metadata_free( vcs_metadata );
    mwalib_voltage_context_free( vcs_context );

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


