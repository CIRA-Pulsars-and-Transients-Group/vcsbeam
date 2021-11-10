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
    bool  empty_lines;     // Insert empty lines between channels in output
    int   nchans;          // The number of frequencies to be used in calculation
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
    MetafitsContext  *obs_context  = NULL;
    MetafitsMetadata *obs_metadata = NULL;
    get_mwalib_metafits_metadata( opts.metafits, &obs_metadata, &obs_context );

    // Set the channels to use in calculation
    if (opts.nchans <= 0)
        opts.nchans = obs_metadata->num_metafits_coarse_chans;

    uint32_t freq_hz; // Used for loop iterator
    double BW_hz = (double)(obs_metadata->obs_bandwidth_hz); // Total bandwidth
    double   bw_hz = BW_hz / opts.nchans; // Make floating point to allow fractional-Hz divisions
    double   flo_hz = obs_metadata->centre_freq_hz - BW_hz/2.0;
    double   freq_hz_start = flo_hz + bw_hz/2.0;

    // Get a "beam geometry" struct (and other variables) ready
    struct beam_geom bg, arrf_bg;
    double mjd;
    double az, za; // Shorthand for azimuth and zenith angle
    double IQUV[4];
    double array_factor;

    primary_beam pb;
    cuDoubleComplex *J = NULL; // For the FEE beam
    uintptr_t coarse_chan_idx = 0; // <-- just a dummy for initially setting up the primary beam struct
    uintptr_t npointings = 1;
    create_primary_beam( &pb, obs_metadata, coarse_chan_idx, npointings );

    // This program assumes no dead dipoles
    uint32_t *delays = obs_metadata->delays;
    double amps[obs_metadata->num_delays];
    uintptr_t i;
    for (i = 0; i < obs_metadata->num_delays; i++)
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
            "# 10 Array factor\n"
            "# 11 FEE beam (PX) real\n"
            "# 12 FEE beam (PX) imag\n"
            "# 13 FEE beam (PY) real\n"
            "# 14 FEE beam (PY) imag\n"
            "# 15 FEE beam (QX) real\n"
            "# 16 FEE beam (QX) imag\n"
            "# 17 FEE beam (QY) real\n"
            "# 18 FEE beam (QY) imag\n"
            "#\n",
            obs_metadata->sched_start_mjd,
            obs_metadata->az_deg, obs_metadata->za_deg,
            obs_metadata->ra_tile_pointing_deg, obs_metadata->dec_tile_pointing_deg,
            ra_hours * PAL__DH2R * PAL__DR2D, dec_degs
           );

    // Loop over the coarse channels
    uintptr_t c, t;
    for (c = 0; c < opts.nchans; c++)
    {
        // Get the frequency of this channel
        freq_hz = (uint32_t)round(freq_hz_start + c*bw_hz);

        // Loop over the gps seconds
        for (t = opts.time_stride/2; t < obs_metadata->num_metafits_timesteps; t += opts.time_stride)
        {
            // Calculate the beam geometry for the requested pointing
            mjd = obs_metadata->sched_start_mjd + (double)t/86400.0;
            calc_beam_geom( ra_hours, dec_degs, mjd, &bg );
            az = bg.az;
            za = PAL__DPIBY2 - bg.el;

            array_factor = 1.0;
            if (opts.do_array_factor)
            {
                calc_beam_geom( arrf_ra_hours, arrf_dec_degs, mjd, &arrf_bg );
                array_factor = calc_array_factor( obs_metadata, freq_hz, &bg, &arrf_bg );
            }

            calc_normalised_beam_response( pb.beam, az, za, freq_hz, delays, amps, IQUV, &J );

            // Print out the results
            fprintf( opts.fout, "%lu %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    t, freq_hz/1e6,
                    az*PAL__DR2D, za*PAL__DR2D,
                    palDsep( az, bg.el,
                        obs_metadata->az_deg*PAL__DD2R, obs_metadata->alt_deg*PAL__DD2R ) * PAL__DR2D,
                    IQUV[0],
                    IQUV[1],
                    IQUV[2],
                    IQUV[3],
                    array_factor,
                    cuCreal( J[0] ), cuCimag( J[0] ),
                    cuCreal( J[1] ), cuCimag( J[1] ),
                    cuCreal( J[2] ), cuCimag( J[2] ),
                    cuCreal( J[3] ), cuCimag( J[3] )
                   );

            free( J );
        }

        // Insert a blank line in the output, to delimit different frequencies
        if (opts.empty_lines)
            fprintf( opts.fout, "\n" );
    }

    // Close output file, if necessary
    if (opts.fout != stderr)
    {
        fclose( opts.fout );
    }

    free_primary_beam( &pb );

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
            "\t-e, --empty-lines          Insert empty lines between channels in output [default: off]\n"
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
    opts->empty_lines     = false;
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
                {"RA",              required_argument, 0, 'r'},
                {"RA-tied",         required_argument, 0, 'R'},
                {"time-stride",     required_argument, 0, 't'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv, "c:d:D:ehm:o:r:R:t:", long_options, &option_index);

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
                case 'r':
                    opts->ra_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->ra_str, optarg );
                    break;
                case 'R':
                    opts->arrf_ra_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->arrf_ra_str, optarg );
                    break;
                case 't':
                    opts->time_stride = atoi(optarg);
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

