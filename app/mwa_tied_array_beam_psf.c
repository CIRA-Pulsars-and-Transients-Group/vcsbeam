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
    double height_deg;     // Height of image in degrees
    double width_deg;      // Width of image in degrees
    int height_pixels;     // Height of image in pixels
    int width_pixels;      // Width of image in pixels
    double freq_MHz;       // The frequency in MHz
    FILE *fout;            // Where to put the output (default STDOUT)
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

    // Set up the grid coordinates
    double width_hours = opts.width_deg / cos( dec_degs*D2R ) * D2R*R2H;
    double dX = width_hours / opts.width_pixels;
    double dY = opts.height_deg / opts.height_pixels;
    double X0 = ra_hours - width_hours/2.0 + dX/2.0;
    double Y0 = dec_degs - opts.height_deg/2.0 + dY/2.0;

    // Get the metadata for the selected observation
    MetafitsContext  *obs_context  = NULL;
    MetafitsMetadata *obs_metadata = NULL;
    get_mwalib_metafits_metadata( opts.metafits, &obs_metadata, &obs_context );

    // If no frequency was selected, get the obs centre frequency
    uint32_t freq_hz;
    if (opts.freq_MHz == 0.0)
        freq_hz = obs_metadata->centre_freq_hz;
    else
        freq_hz = (uint32_t)(opts.freq_MHz * 1e6);

    // Get a "beam geometry" struct (and other variables) ready
    struct beam_geom bg, arrf_bg;
    double mjd;
    double az, za; // Shorthand for azimuth and zenith angle
    double IQUV[4];
    double array_factor;

    primary_beam pb;
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
            "# Point spread function for an MWA tied-array beam\n"
            "# ------------------------------------------------\n"
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
            "# Frequency:            %f MHz\n"
            "#\n"
            "# 1            2            3\n"
            "# RA (hours) | Dec (degs) | Array_factor\n",
            obs_metadata->sched_start_mjd,
            obs_metadata->az_deg, obs_metadata->za_deg,
            obs_metadata->ra_tile_pointing_deg, obs_metadata->dec_tile_pointing_deg,
            ra_hours * H2R * R2D, dec_degs,
            opts.freq_MHz
           );

    // Calculate the beam geometry for the requested pointing
    mjd = obs_metadata->sched_start_mjd;
    calc_beam_geom( ra_hours, dec_degs, mjd, &bg );
    az = bg.az;
    za = PIBY2 - bg.el;

    // Loop over RA
    int X_idx, Y_idx;
    double X, Y;
    cuDoubleComplex *J=  NULL;
    for (X_idx = 0; X_idx < opts.width_pixels; X_idx++)
    {
        X = X0 + X_idx*dX;

        // Loop over Dec
        for (Y_idx = 0; Y_idx < opts.height_pixels; Y_idx++)
        {
            Y = Y0 + Y_idx*dY;

            calc_beam_geom( X, Y, mjd, &arrf_bg );
            array_factor = calc_array_factor( obs_metadata, freq_hz, &arrf_bg, &bg );

            calc_normalised_beam_response( pb.beam, az, za, freq_hz, delays, amps, IQUV, &J, true );

            // Print out the results
            fprintf( opts.fout, "%f %f %f\n",
                    X, Y, array_factor*IQUV[0] );

            free( J );
        }

        // Insert a blank line in the output, to delimit different rows
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
            "\t-f, --freq_MHz=MHZ         Frequency in MHz [default: mid-frequency of observation]\n"
            "\t-x, --width=DEG            Width of image in degrees [default: 1.0]\n"
            "\t-X, --npixelsx=PIXELS      Width of image in pixels [default: 100]\n"
            "\t-y, --height=DEG           Height of image in degrees [default: 1.0]\n"
            "\t-Y, --npixelsy=PIXELS      Height of image in pixels [default: 100]\n"
            "\t-o, --outfile=FILENAME     Write the results to FILENAME [default is to write to STDOUT]\n"
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
    opts->height_deg      = 1.0;     // Height of image in degrees
    opts->width_deg       = 1.0;     // Width of image in degrees
    opts->height_pixels   = 100;     // Height of image in pixels
    opts->width_pixels    = 100;     // Width of image in pixels
    opts->freq_MHz        = 0.0;     // The frequency in MHz (0.0 means get obs centre frequency)
    opts->fout            = stdout;

    if (argc > 1)
    {
        int c;
        while (1)
        {
            static struct option long_options[] = {
                {"metafits",        required_argument, 0, 'm'},
                {"RA",              required_argument, 0, 'r'},
                {"Dec",             required_argument, 0, 'd'},
                {"width",           required_argument, 0, 'x'},
                {"npixelsx",        required_argument, 0, 'X'},
                {"height",          required_argument, 0, 'y'},
                {"npixelsy",        required_argument, 0, 'Y'},
                {"freq_MHz",        required_argument, 0, 'f'},
                {"outfile",         required_argument, 0, 'o'},
                {"help",            required_argument, 0, 'h'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv, "d:f:hm:o:r:x:X:y:Y:", long_options, &option_index);

            if (c == -1)
                break;

            switch( c )
            {
                case 'd':
                    opts->dec_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->dec_str, optarg );
                    break;
                case 'f':
                    opts->freq_MHz = atof( optarg );
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
                case 'x':
                    opts->width_deg = atof( optarg );
                    break;
                case 'y':
                    opts->height_deg = atof( optarg );
                    break;
                case 'X':
                    opts->width_pixels = atoi( optarg );
                    break;
                case 'Y':
                    opts->height_pixels = atoi( optarg );
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
    // and that all the numbers are sensible
    assert( opts->metafits       != NULL );
    assert( opts->ra_str         != NULL );
    assert( opts->dec_str        != NULL );
    assert( opts->width_deg      > 0.0 );
    assert( opts->height_deg     > 0.0 );
    assert( opts->width_pixels   > 0   );
    assert( opts->height_pixels  > 0   );

}

