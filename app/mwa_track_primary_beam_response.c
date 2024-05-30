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

// Non-standard dependencies
#include <cuComplex.h>
#include <mwalib.h>
#include <star/pal.h>
#include <star/palmac.h>
#include <mwa_hyperbeam.h>

// Local includes
#include "vcsbeam.h"

#define NCOMPLEXELEMENTS 4

const double Isky[] = { 0.5, 0. , 0. ,  0. , 0. , 0. ,  0.5, 0. };
const double Qsky[] = { 0.5, 0. , 0. ,  0. , 0. , 0. , -0.5, 0. };
const double Usky[] = { 0. , 0. , 0.5,  0. , 0.5, 0. ,  0. , 0. };
const double Vsky[] = { 0. , 0. , 0. , -0.5, 0  , 0.5,  0. , 0. };
const double *sky[] = { Isky, Qsky, Usky, Vsky };

struct mwa_track_primary_beam_response_opts {
    char *metafits;        // filename of the metafits file
    char *ra_str;          // String representing the RA
    char *dec_str;         // String representing the Dec
    char *arrf_ra_str;     // String representing the RA of a tied-array beam
    char *arrf_dec_str;    // String representing the Dec of a tied-array beam
    FILE *fout;            // Where to put the output (default STDOUT)
    bool  do_array_factor; // Whether to calculate the array factor or not
};

void usage();
void mwa_track_primary_beam_response_parse_cmdline(
        int argc, char **argv, struct mwa_track_primary_beam_response_opts *opts );
void calc_normalised_beam_response( FEEBeam *beam, double az, double za, double freq_hz, uint32_t *delays, double *amps, double *IQUV, cuDoubleComplex *J );

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

    // Get a "beam geometry" struct (and other variables) ready
    struct beam_geom bg, arrf_bg;
    double mjd;
    double az, za; // Shorthand for azimuth and zenith angle
    uint32_t freq_hz;
    double IQUV[4];
    double array_factor;

    primary_beam pb;
    cuDoubleComplex *J = malloc( 4*sizeof(cuDoubleComplex) );
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
            "# 1        2                3            4                5                             6  7  8  9  10\n"
            "# Seconds  Frequency_(MHz)  Azimuth_(°)  ZenithAngle_(°)  Distance_from_tile_centre_(°) I  Q  U  V  Array_factor\n",
            obs_metadata->sched_start_mjd,
            obs_metadata->az_deg, obs_metadata->za_deg,
            obs_metadata->ra_tile_pointing_deg, obs_metadata->dec_tile_pointing_deg,
            ra_hours * PAL__DH2R * PAL__DR2D, dec_degs
           );

    // Loop over the coarse channels
    uintptr_t c, t;
    for (c = 0; c < obs_metadata->num_metafits_coarse_chans; c++)
    {
        // Get the frequency of this channel
        freq_hz = obs_metadata->metafits_coarse_chans[c].chan_centre_hz;

        // Loop over the gps seconds
        for (t = 0; t < obs_metadata->num_metafits_timesteps; t++)
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

            calc_normalised_beam_response( pb.beam, az, za, freq_hz, delays, amps, IQUV, J );

            // Print out the results
            fprintf( opts.fout, "%lu %f %f %f %f %f %f %f %f %f\n",
                    t, freq_hz/1e6,
                    az*PAL__DR2D, za*PAL__DR2D,
                    palDsep( az, bg.el,
                        obs_metadata->az_deg*PAL__DD2R, obs_metadata->alt_deg*PAL__DD2R ) * PAL__DR2D,
                    IQUV[0],
                    IQUV[1],
                    IQUV[2],
                    IQUV[3],
                    array_factor );
        }

        // Insert a blank line in the output, to delimit different frequencies
        fprintf( opts.fout, "\n" );
    }

    free( J );

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

    if (argc > 1)
    {
        int c;
        while (1)
        {
            static struct option long_options[] = {
                {"metafits",        required_argument, 0, 'm'},
                {"RA",              required_argument, 0, 'r'},
                {"Dec",             required_argument, 0, 'd'},
                {"RA-tied",         required_argument, 0, 'R'},
                {"Dec-tied",        required_argument, 0, 'D'},
                {"outfile",         required_argument, 0, 'o'},
                {"help",            required_argument, 0, 'h'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv, "d:D:hm:o:r:R:", long_options, &option_index);

            if (c == -1)
                break;

            switch( c )
            {
                case 'd':
                    opts->dec_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->dec_str, optarg );
                    break;
                case 'D':
                    opts->arrf_dec_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->arrf_dec_str, optarg );
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

    if (opts->arrf_ra_str && opts->arrf_dec_str)
    {
        opts->do_array_factor = true;
    }
}

void calc_normalised_beam_response( FEEBeam *beam, double az, double za, double freq_hz, uint32_t *delays, double *amps, double *IQUV, cuDoubleComplex *J )
{
    uint8_t iauOrder = 0; //Boolean 0= don't set Jones matrix to be iau order but instead mwa order in CalcJones

    uint32_t numAmps = 16; //number of dipole gains used (16 or 32)

    int32_t errInt = 0; //exit code integer for calcJones

    double *arrayLatitudeRad = NULL;

    cuDoubleComplex JH[NCOMPLEXELEMENTS];
    cuDoubleComplex sky_x_JH[NCOMPLEXELEMENTS];
    cuDoubleComplex coherency[NCOMPLEXELEMENTS];
    double P[NCOMPLEXELEMENTS]; // (Real-valued) parallactic angle correction matrix

    // Calculate the parallactic angle correction
    parallactic_angle_correction( P, MWA_LATITUDE_RADIANS, az, za );

    // Calculate the primary beam for this channel, in this direction
    int zenith_norm = 1;
    errInt = calc_jones(beam, az, za, freq_hz, delays, amps, numAmps, zenith_norm, arrayLatitudeRad, iauOrder, (double *)J);

    if (errInt !=0)
    {             
        printf("Hyperbeam failed with exit code %d\n", EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    mult2x2d_RxC( P, J, J );

    // Convert this jones matrix to Stokes parameters for I, Q, U, V skies
    int stokes;
    for (stokes = 0; stokes < 4; stokes++)
    {
        calc_hermitian( J, JH );
        mult2x2d( (cuDoubleComplex *)sky[stokes], JH, sky_x_JH );
        mult2x2d( J, sky_x_JH, coherency );

        if (stokes == 0) // Stokes I
            IQUV[0] = 0.5*cuCreal( cuCadd( coherency[0], coherency[3] ) );
        else if (stokes == 1) // Stokes Q
            IQUV[1] = 0.5*cuCreal( cuCsub( coherency[0], coherency[3] ) );
        else if (stokes == 2) // Stokes U
            IQUV[2] =  cuCreal( coherency[1] );
        else  // if (stokes == 3) // Stokes V
            IQUV[3] = -cuCimag( coherency[1] );
    }

    free( J );
}
