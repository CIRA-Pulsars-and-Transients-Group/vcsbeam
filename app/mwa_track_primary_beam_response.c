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
#include <time.h>

// Local includes
#include "vcsbeam.h"

struct mwa_track_primary_beam_response_opts {
    char *metafits;        // filename of the metafits file
    char *ra_str;          // String representing the RA
    char *dec_str;         // String representing the Dec
    char *arrf_ra_str;     // String representing the RA of a tied-array beam
    char *arrf_dec_str;    // String representing the Dec of a tied-array beam
    FILE *fout;            // Where to put the output (default STDOUT)
    bool  do_array_factor; // Whether to calculate the array factor or not
    int   time_stride;     // Output one measurement per TIME_STRIDE seconds
    int   start_time;      // Only go up to START_TIME seconds
    int   end_time;        // Only go up to (but not including) END_TIME seconds
    bool  empty_lines;     // Insert empty lines between channels in output
    int   nchans;          // The number of frequencies to be used in calculation
    bool  apply_pa_correction; // Apply the parallactic angle correction
};

void usage();
void mwa_track_primary_beam_response_parse_cmdline(
        int argc, char **argv, struct mwa_track_primary_beam_response_opts *opts );

int main(int argc, char **argv)
{
    // Parse command line arguments
    struct mwa_track_primary_beam_response_opts opts;
    mwa_track_primary_beam_response_parse_cmdline( argc, argv, &opts );

    // Parse the RA and Dec
    double ra_hours = parse_ra( opts.ra_str );
    double dec_degs = parse_dec( opts.dec_str );

    double arrf_ra_hours = parse_ra( opts.arrf_ra_str );
    double arrf_dec_degs = parse_dec( opts.arrf_dec_str );

    // Get the metadata for the selected observation
    vcsbeam_context vm;
    vmLoadObsMetafits( &vm, opts.metafits );

    // Set the channels to use in calculation
    if (opts.nchans <= 0)
        opts.nchans = vm.obs_metadata->num_metafits_coarse_chans;

    // Set start and end times, if not explicitly given on command line, or if
    // invalid start/end times are given
    if (opts.start_time < 0)              opts.start_time = 0;
    if (opts.end_time < opts.start_time)  opts.end_time   = vm.obs_metadata->num_metafits_timesteps;

    uint32_t freq_hz; // Used for loop iterator
    double BW_hz = (double)(vm.obs_metadata->obs_bandwidth_hz); // Total bandwidth
    double   bw_hz = BW_hz / opts.nchans; // Make floating point to allow fractional-Hz divisions
    double   flo_hz = vm.obs_metadata->centre_freq_hz - BW_hz/2.0;
    double   freq_hz_start = flo_hz + bw_hz/2.0;

    // Get a "beam geometry" struct (and other variables) ready
    beam_geom bg, arrf_bg;
    double mjd;
    double az, za; // Shorthand for azimuth and zenith angle
    double IQUV[4];
    double array_factor;

    char coord1[8], coord2[8];
    if (opts.apply_pa_correction)
    {
        sprintf( coord1, "Y" );
        sprintf( coord2, "X" );
    }
    else
    {
        sprintf( coord1, "θ" );
        sprintf( coord2, "ϕ" );
    }

    gpuDoubleComplex *J = malloc( 4*sizeof(gpuDoubleComplex) ); // For the FEE beam
    vm.npointing = 1;
    vm.coarse_chan_idx = 0; // <-- just a dummy for initially setting up the primary beam struct
    vmCreatePrimaryBeam( &vm );

    // This program assumes no dead dipoles
    uint32_t *delays = vm.obs_metadata->delays;
    double amps[vm.obs_metadata->num_delays];
    uintptr_t i;
    for (i = 0; i < vm.obs_metadata->num_delays; i++)
        amps[i] = 1.0;

    // Write out the output header
    fprintf( opts.fout,
            "# Zenith-normalised Stokes response for given primary beam and tied-array beam\n"
            "# ----------------------------------------------------------------------------\n"
            "# vcsbeam %s\n"
            "#    ",
            VCSBEAM_VERSION );
    for (i = 0; i < (uintptr_t)argc; i++)
        fprintf( opts.fout, " %s", argv[i] );
    fprintf( opts.fout,
            "\n"
            "# MJD:                  %f\n"
            "# Tile pointing centre: Azimuth = %f°; Zenith angle = %f°\n"
            "# (Initial) tile pos:   RA      = %f°; Dec = %f°\n"
            "# Tied array position:  RA      = %f°; Dec = %f°\n"
            "#\n"
            "# Columns:\n"
            "# 1  Seconds into observation\n"
            "# 2  Frequency                 (MHz)\n"
            "# 3  Azimuth                   (°)\n"
            "# 4  Zenith Angle              (°)\n"
            "# 5  Distance from tile centre (°)\n"
            "# 6  I                         (Power, a.u.)\n"
            "# 7  Q                         (Power, a.u.)\n"
            "# 8  U                         (Power, a.u.)\n"
            "# 9  V                         (Power, a.u.)\n"
            "# 10 Array factor\n",
            vm.obs_metadata->sched_start_mjd,
            vm.obs_metadata->az_deg, vm.obs_metadata->za_deg,
            vm.obs_metadata->ra_tile_pointing_deg, vm.obs_metadata->dec_tile_pointing_deg,
            ra_hours * PAL__DH2R * PAL__DR2D, dec_degs
           );

    for (i = 0; i < 8; i++)
        fprintf( opts.fout,
            "# %d FEE beam (%c%s) %s\n",
            i+11, (i/4 ? 'P' : 'Q'), (i%4 < 2 ? coord1 : coord2), (i%2 ? "real" : "imag") );

    fprintf( opts.fout, "#\n" );

    // Loop over the coarse channels
    uintptr_t c, t;
    for (c = 0; c < opts.nchans; c++)
    {
        // Get the frequency of this channel
        freq_hz = (uint32_t)round(freq_hz_start + c*bw_hz);

        // Loop over the gps seconds
        for (t = opts.start_time; t < opts.end_time; t += opts.time_stride)
        {
            // Calculate the beam geometry for the requested pointing
            mjd = vm.obs_metadata->sched_start_mjd + (double)t/86400.0;
            calc_beam_geom( ra_hours, dec_degs, mjd, &bg );
            az = bg.az;
            za = PAL__DPIBY2 - bg.el;

            array_factor = 1.0;
            if (opts.do_array_factor)
            {
                calc_beam_geom( arrf_ra_hours, arrf_dec_degs, mjd, &arrf_bg );
                array_factor = calc_array_factor( vm.obs_metadata, freq_hz, &bg, &arrf_bg );
            }

            calc_normalised_beam_response( vm.pb.beam, az, za, freq_hz, delays, amps, IQUV, J, opts.apply_pa_correction );

            // Print out the results
            fprintf( opts.fout, "%lu %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    t, freq_hz/1e6,
                    az*PAL__DR2D, za*PAL__DR2D,
                    palDsep( az, bg.el,
                        vm.obs_metadata->az_deg*PAL__DD2R, vm.obs_metadata->alt_deg*PAL__DD2R ) * PAL__DR2D,
                    IQUV[0],
                    IQUV[1],
                    IQUV[2],
                    IQUV[3],
                    array_factor,
                    gpuCreal( J[0] ), gpuCimag( J[0] ),
                    gpuCreal( J[1] ), gpuCimag( J[1] ),
                    gpuCreal( J[2] ), gpuCimag( J[2] ),
                    gpuCreal( J[3] ), gpuCimag( J[3] )
                   );

        }

        // Insert a blank line in the output, to delimit different frequencies
        if (opts.empty_lines)
            fprintf( opts.fout, "\n" );
    }

    free( J );

    // Close output file, if necessary
    if (opts.fout != stderr)
    {
        fclose( opts.fout );
    }

    free_primary_beam( &vm.pb );

    // Exit gracefully
    return EXIT_SUCCESS;
}

void usage()
{
    printf( "\nusage: mwa_track_primary_beam_response [OPTIONS]\n");

    printf( "\nOPTIONS\n\n"
            "\t-m, --metafits=FILE        FILE is the metafits file for the target observation (required)\n"
            "\t-r, --RA=HH:MM:SS          Tied-array pointing direction, right ascension (required)\n"
            "\t-d, --Dec=DD:MM:SS         Tied-array pointing direction, declination (required)\n"
            "\t-R, --RA-tied=HH:MM:SS     Direction for array factor calculation, right ascension\n"
            "\t-D, --Dec-tied=DD:MM:SS    Direction for array factor calculation, declination\n"
            "\t-c, --num-chans=NCHANS     The number of frequency channels to use in calculation [default: number of coarse channels in obs]\n"
            "\t-t, --time-stride=NSECONDS Output one measurement every NSECONDS seconds [default: 1]\n"
            "\t-S, --start-time=NSECONDS  Start at NSECONDS (from start of observation) [default: start at beginning of observation]\n"
            "\t-T, --end-time=NSECONDS    Only go up to (but not including) NSECONDS (from start of observation) "
                                         "[default: go till end of observation]\n"
            "\t-e, --empty-lines          Insert empty lines between channels in output [default: off]\n"
            "\t-P, --apply-pa             Apply the parallactic angle correction [default: off]\n"
            "\t-o, --outfile=FILENAME     Write the results to FILENAME [default is to write to STDOUT\n"
            "\t-h, --help                 Write this help message and exit\n\n"
          );
}

void mwa_track_primary_beam_response_parse_cmdline(
        int argc, char **argv, struct mwa_track_primary_beam_response_opts *opts )
{
    // Set defaults
    opts->metafits        = NULL;
    opts->ra_str          = NULL;
    opts->dec_str         = NULL;
    opts->arrf_ra_str     = NULL;
    opts->arrf_dec_str    = NULL;
    opts->fout            = stdout;
    opts->do_array_factor = false;
    opts->time_stride     = 1;
    opts->start_time      = -1;
    opts->end_time        = -1;
    opts->empty_lines     = false;
    opts->apply_pa_correction = false;
    opts->nchans          = -1; // "Default" value to indicate (later) to use number of coarse channels

    if (argc > 1)
    {
        int c;
        while (1)
        {
            static struct option long_options[] = {
                {"num-chans",       required_argument, 0, 'c'},
                {"Dec",             required_argument, 0, 'd'},
                {"Dec-tied",        required_argument, 0, 'D'},
                {"empty-lines",     no_argument,       0, 'e'},
                {"help",            no_argument,       0, 'h'},
                {"metafits",        required_argument, 0, 'm'},
                {"outfile",         required_argument, 0, 'o'},
                {"apply-pa",        no_argument,       0, 'P'},
                {"RA",              required_argument, 0, 'r'},
                {"RA-tied",         required_argument, 0, 'R'},
                {"time-stride",     required_argument, 0, 't'},
                {"start-time",      required_argument, 0, 'S'},
                {"end-time",        required_argument, 0, 'T'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv, "c:d:D:ehm:o:Pr:R:S:t:T:", long_options, &option_index);

            if (c == -1)
                break;

            switch( c )
            {
                case 'c':
                    opts->nchans = atoi(optarg);
                    break;
                case 'd':
                    opts->dec_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->dec_str, optarg );
                    break;
                case 'D':
                    opts->arrf_dec_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->arrf_dec_str, optarg );
                    break;
                case 'e':
                    opts->empty_lines = true;
                    break;
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'm':
                    opts->metafits = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->metafits, optarg );
                    break;
                case 'o':
                    opts->fout = fopen( optarg, "w" );
                    if (opts->fout == NULL)
                    {
                        fprintf( stderr, "error: could not open file '%s' for writing\n", optarg );
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'P':
                    opts->apply_pa_correction = true;
                    break;
                case 'r':
                    opts->ra_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->ra_str, optarg );
                    break;
                case 'R':
                    opts->arrf_ra_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->arrf_ra_str, optarg );
                    break;
                case 'S':
                    opts->start_time = atoi(optarg);
                    break;
                case 't':
                    opts->time_stride = atoi(optarg);
                    break;
                case 'T':
                    opts->end_time = atoi(optarg);
                    break;
                default:
                    fprintf( stderr, "error: unrecognised option '%s'\n", optarg );
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

    // Check that all the required arguments have been provided
    assert( opts->metafits       != NULL );
    assert( opts->ra_str         != NULL );
    assert( opts->dec_str        != NULL );
    assert( opts->time_stride     > 0    );

    if (opts->arrf_ra_str && opts->arrf_dec_str)
    {
        opts->do_array_factor = true;
    }
}

