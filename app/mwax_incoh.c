// Standard library
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
//#include <mpi.h>

#include <mwalib.h>

#define ERROR_MESSAGE_LEN  1024
#define PROCESS_UNTIL_END  -1

#define VSPVB  64000  /* "Vcsmwax Samples Per Voltage Block"
                         (i.e. the size of the TIME dimension in a voltage block */
#define vMWAX_IDX(s,i,ni) (((s/VSPVB)*(ni) + (i))*VSPVB + (s%VSPVB))

struct cmd_line_opts {
    int skip;                // Number of seconds to skip at the beginning
    int duration;            // How many seconds to process
    char *datadir;           // The path to where the recombined data live
    char *metafits;          // Path to the metafits file
    char *outfile;           // Base name of the output PSRFITS file
    char *flags;             // Tiles to be flagged
};

void make_incoh_beam_parse_cmdline(
        int argc, char **argv, struct cmd_line_opts *opts );
void usage(FILE *f);

int main(int argc, char *argv[])
{
    struct cmd_line_opts opts;
    make_incoh_beam_parse_cmdline( argc, argv, &opts );

    char error_message[ERROR_MESSAGE_LEN];

    MetafitsContext  *obs_context;
    MetafitsMetadata *obs_metadata;
    VoltageContext   *vcs_context;
    VoltageMetadata  *vcs_metadata;

    int coarse_chan_idx = 3; // i.e. rcvr chan 160 (the 0th channel is 157)

    // mwalib contexts
    if (mwalib_metafits_context_new2(opts.metafits, &obs_context, error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    if (mwalib_metafits_metadata_get(obs_context, NULL, NULL, &obs_metadata, error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    TimeStep timestep;
    for (int i = 0; i < obs_metadata->num_metafits_timesteps; i++)
    {
        timestep = obs_metadata->metafits_timesteps[i];
        printf("%lu\n", timestep.gps_time_ms);
    }

    // Use metafits info to work out start and end times, and number of timesteps
    int first_gps_second = obs_metadata->sched_start_gps_time_ms/1000 + opts.skip;
    // Calculate the first good GPS second to start on
    // mwalib gives the first good GPS millisecond, so for GPS second 1111111111,
    //    GPS millisecond  ->  GPS second
    //     1111111110998   ->  1111111111
    //     1111111110999   ->  1111111111
    //     1111111111000   ->  1111111111
    //     1111111111001   ->  1111111112
    //     1111111111002   ->  1111111112
    // In integer arithmetic, this is achieved via
    //          ms         ->  (ms - 1) / 1000 + 1
    int first_good_gps_second = (obs_metadata->good_time_gps_ms - 1) / 1000 + 1;
    // Calculate which GPS second to stop at
    // mwalib gives the secheduled end in GPS ms, so the mapping goes like
    //    GPS millisecond  ->  Have to stop before this GPS second
    //     1111111110998   ->  1111111110
    //     1111111110999   ->  1111111110
    //     1111111111000   ->  1111111111
    //     1111111111001   ->  1111111111
    int stop_gps_second = obs_metadata->sched_end_gps_time_ms / 1000;
    //if last_gps_second = 
    // Get filename:
    const char *voltage_file = malloc(1024);

    if (mwalib_metafits_get_expected_volt_filename(
                obs_context,
                0,
                coarse_chan_idx,
                voltage_file,
                1024,
                error_message,
                ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot get expected volt filename: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    if (mwalib_voltage_context_new(opts.metafits, &voltage_file, 1, &vcs_context, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    if (mwalib_voltage_metadata_get(vcs_context, &vcs_metadata, error_message, ERROR_MESSAGE_LEN) != EXIT_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot get metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Replace the existing OBS_METADATA with a new one that knows this is a VCS observation
    mwalib_metafits_metadata_free( obs_metadata );
    if (mwalib_metafits_metadata_get( NULL, NULL, vcs_context, &obs_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata from voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create buffer for 1 second's worth of data
    size_t buffer_size = vcs_metadata->num_voltage_blocks_per_second * vcs_metadata->voltage_block_size_bytes;
    char *buffer = malloc(buffer_size);

    int gps_second;
    for (gps_second = first_good_gps_second; gps_second <= 1367261002; gps_second++)
    {
        // Print out read function inputs
        //mwalib_voltage_context_display(vcs_context, error_message, ERROR_MESSAGE_LEN);

        if (mwalib_voltage_context_read_second(vcs_context, gps_second, 1, coarse_chan_idx, buffer, buffer_size, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
        {
            fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s", error_message );
            exit(EXIT_FAILURE);
        }
    }

    // Free memory
    free(buffer);

    return EXIT_SUCCESS;
}

void usage(FILE *f)
{
    fprintf(f, "Incoherent beamformer for MWAX data\n");
    fprintf(f, "===================================\n");
    fprintf(f, "usage: mwax_incoh [options] METAFITS\n\n");
    fprintf(f, "Required arguments:\n");
    fprintf(f, "  METAFITS       The path to the observation's metafits file\n");
    fprintf(f, "Options:\n");
    fprintf(f, "  -d DATADIR     The directory containing the data files. If not given, assumed to be\n");
    fprintf(f, "                 the same directory as the metafits file.\n");
    fprintf(f, "  -f TILENAME[,TILENAME[,...]]\n");
    fprintf(f, "                 A comma-separated list of tilenames to be flagged (i.e. excluded from\n");
    fprintf(f, "                 the incoherent beam)\n");
    fprintf(f, "  -h             Print this help and exit\n");
    fprintf(f, "  -o OUTFILE     Basename of the output PSRFITS file\n");
    fprintf(f, "  -s SECONDS     Skip the first SECONDS seconds of data. If not supplied, start at the \n");
    fprintf(f, "  -t SECONDS     Process SECONDS seconds of data. If not supplied, process until the end.\n");
}

void make_incoh_beam_parse_cmdline(
        int argc, char **argv, struct cmd_line_opts *opts )
{
    // Set defaults for command line options
    opts->skip        = 0;
    opts->duration    = PROCESS_UNTIL_END;
    opts->datadir     = NULL;
    opts->outfile     = NULL;
    opts->metafits    = NULL;
    opts->flags       = NULL;

    int c, opt;
    while ((opt = getopt(argc, argv, "d:f:ho:s:t:")) != -1)
    {
        switch(opt) {

            case 'd':
                opts->datadir = strdup(optarg); // Duplicate because this will be freed regardless
                break;
            case 'f':
                opts->flags = optarg;
                break;
            case 'h':
                usage(stdout);
                exit(EXIT_SUCCESS);
                break;
            case 'o':
                opts->outfile = optarg;
                break;
            case 's':
                opts->skip = atoi(optarg);
                break;
            case 't':
                opts->duration = atoi(optarg);
                break;
            default:
                fprintf(stderr, "unrecognised option '%s'\n", optarg);
                usage(stderr);
                exit(EXIT_FAILURE);
        }
    }

    if (optind >= argc)
    {
        fprintf(stderr, "no metafits supplied\n");
        usage(stderr);
        exit(EXIT_FAILURE);
    }

    opts->metafits = argv[optind];

    // If not supplied, set the default data directory to equal the location
    // of the metafits file
    if (opts->datadir == NULL)
    {
        char *path_end = strrchr(opts->metafits, '/');
        if (path_end != NULL)
        {
            int len = path_end - opts->metafits + 2; // The length of the string to copy, +2 for '/' and '\0'
            opts->datadir = malloc(len);
            strncpy(opts->datadir, opts->metafits, len);
            opts->datadir[len-1] = '\0';
        }
    }
}
