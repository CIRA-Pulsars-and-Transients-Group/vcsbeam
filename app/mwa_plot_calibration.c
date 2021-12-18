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

// Non-standard dependencies
#include <mwalib.h>

// Local includes
#include "vcsbeam.h"

struct make_plot_calibrate_opts {
    char      *metafits; // filename of the metafits file
    uintptr_t  ncoarse_chans;
    char      *custom_flags;
    char      *zero_string;

    // Calibration options
    char              *cal_metafits;     // Filename of the metafits file
    char              *caldir;           // Location of calibration data
    int                cal_type;         // Either RTS or OFFRINGA
    char              *ref_ant;          // Reference antenna for calibration phases
    double             phase_offset;     // Rotate the phase of Y by m*freq + c, where
    double             phase_slope;      //   m = phase_slope (rad/Hz)
                                         //   c = phase_offset (rad)
    bool               custom_pq_correction; // Set to true if phase_offset and phase_slope are to be applied
    bool               keep_cross_terms; // Include PQ and QP of calibration Jones matrices
    bool               use_bandpass;     // Use the Bandpass solutions
};

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/

void usage();
void make_plot_calibrate_parse_cmdline( int argc, char **argv,
        struct make_plot_calibrate_opts *opts );

/********
 * MAIN *
 ********/

int main(int argc, char **argv)
{
    // Parse command line arguments
    struct make_plot_calibrate_opts opts;
    make_plot_calibrate_parse_cmdline( argc, argv, &opts );

    calibration cal;  // Variables for calibration settings
    init_calibration( &cal );

    cal.metafits     = strdup( opts.cal_metafits );
    cal.caldir       = strdup( opts.caldir );
    cal.cal_type     = opts.cal_type;
    cal.ref_ant      = strdup( opts.ref_ant );
    cal.phase_offset = opts.phase_offset;
    cal.phase_slope  = opts.phase_slope;
    cal.custom_pq_correction = opts.custom_pq_correction;
    cal.keep_cross_terms     = opts.keep_cross_terms;
    cal.use_bandpass         = opts.use_bandpass;

    int i; // Generic counter

    // Start a logger for output messages and time-keeping
    logger *log = create_logger( stderr, 0 );

    // Get metadata for obs and cal
    vcsbeam_context vm;
    vmLoadMetafits( opts.metafits, &vm.obs_metadata, &vm.obs_context );
    vmLoadMetafits( cal.metafits, &vm.cal_metadata, &vm.cal_context );

    // Create some "shorthand" variables for code brevity
    uintptr_t nants          = vm.obs_metadata->num_ants;
    uintptr_t nchans         = vm.obs_metadata->num_volt_fine_chans_per_coarse;
    uintptr_t npols          = vm.obs_metadata->num_ant_pols;   // (X,Y)

    // Now, do the following for each coarse channel
    uintptr_t ncoarse_chans = opts.ncoarse_chans;
    if (ncoarse_chans <= 0 || ncoarse_chans > vm.obs_metadata->num_metafits_coarse_chans)
        ncoarse_chans = vm.obs_metadata->num_metafits_coarse_chans;
    uintptr_t Ch, ch; // (Coarse channel idx, fine channel idx)
    cuDoubleComplex *D[ncoarse_chans]; // See Eqs. (27) to (29) in Ord et al. (2019)
    logger *plog = log;

    vmMallocDHost( &vm );

    for (Ch = 0; Ch < ncoarse_chans; Ch++)
    {
        // Read in the calibration solution
        if (cal.cal_type == CAL_RTS)
        {
            vmLoadRTSSolution( &vm, cal.use_bandpass, cal.caldir, Ch );
        }
        else if (cal.cal_type == CAL_OFFRINGA)
        {
            vmLoadOffringaSolution( &vm, Ch, cal.caldir );
        }

        // Copy the solution into the "D" arrays
        D[Ch] = (cuDoubleComplex *)malloc( vm.D_size_bytes );
        memcpy( D[Ch], vm.D, vm.D_size_bytes );

        // Flag antennas that need flagging
        vmSetCustomTileFlags( &vm, opts.custom_flags, &cal );

        // Apply any calibration corrections
        parse_calibration_correction_file( vm.cal_metadata->obs_id, &cal );
        apply_calibration_corrections( &cal, D[Ch], vm.obs_metadata, Ch, plog );

        // After the first coarse channel, don't print out any more messages
        plog = NULL;
    }

    vmFreeDHost( &vm );

    // Print out results
    FILE *fout = stdout;
    fprintf( fout, "# CALIBRATION SOLUTION:\n#\n"
                   "# Created by:\n"
                   "#    " );
    for (i = 0; i < argc; i++)
        fprintf( fout, " %s", argv[i] );
    fprintf( fout, "\n# Column description:\n"
                   "#   Eight columns per antenna, where the eight columns are the magnitudes and phases of\n"
                   "#   the four Jones matrix elements:\n"
                   "#     [ Dqq Dqp ]\n"
                   "#     [ Dpq Dpp ]\n"
                   "# e.g.\n"
                   "#   Tile1_abs(Dqq), Tile1_arg(Dqq), Tile1_abs(Dqp), Tile1_arg(Dqp), ..., Tile2_abs(Dqq), ...\n"
                   "#\n# The ordered tile names are (reading rows first):" );
    uintptr_t ant;
    for (ant = 0; ant < vm.obs_metadata->num_ants; ant++)
    {
        if (ant % 8 == 0)
            fprintf( fout, "\n#  " );
        fprintf( fout, "  %-10s", vm.obs_metadata->antennas[ant].tile_name );
    }
    fprintf( fout, "\n#\n# DATA:\n#\n" );

    int didx; // idx of element in D array
    int pol1, pol2;

    for (Ch = 0; Ch < ncoarse_chans; Ch++)
    {
        for (ch = 0; ch < nchans; ch++)
        {
            for (ant = 0; ant < vm.obs_metadata->num_ants; ant++)
            {
                // Special output if all elements are zero (if user requested)
                didx = J_IDX(ant,ch,0,0,nchans,npols);
                if (opts.zero_string != NULL && is2x2zero( &(D[Ch][didx]) ))
                {
                    for (i = 0; i < npols*npols*2; i++)
                        fprintf( fout, "%s ", opts.zero_string );
                }
                else
                {
                    for (pol1 = 0; pol1 < npols; pol1++)
                        for (pol2 = 0; pol2 < npols; pol2++)
                        {
                            didx = J_IDX(ant,ch,pol1,pol2,nchans,npols);
                            // The jones matrix element to be printed is D[Ch][didx]
                            fprintf( fout, "%e %e ",
                                    cuCabs(D[Ch][didx]),                                // magnitude
                                    atan2( cuCimag(D[Ch][didx]), cuCreal(D[Ch][didx]) ) // phase
                                   );
                        }
                }
            }

            // Print a new line here
            fprintf( fout, "\n" );
        }
    }

    // Free up memory
    for (Ch = 0; Ch < ncoarse_chans; Ch++)
    {
        free( D[Ch] );
    }

    free( cal.caldir           );
    free( opts.metafits        );

    if (opts.custom_flags != NULL)
        free( opts.custom_flags );

    // Free mwalib structs
    mwalib_metafits_context_free( vm.obs_context );
    mwalib_metafits_context_free( vm.cal_context );
    mwalib_metafits_metadata_free( vm.obs_metadata );
    mwalib_metafits_metadata_free( vm.cal_metadata );

    destroy_logger( log );

    return EXIT_SUCCESS;
}


void usage()
{
    printf( "\nusage: mwa_plot_calibration [OPTIONS]\n");

    printf( "\nOPTIONS (RTS)\n\n"
            "\t-m, --metafits=FILE        FILE is the metafits file for the target observation\n"
            "\t-B, --bandpass             Use the Bandpass calibrations (as well as the DIJones solutions) [default: off]\n"
            "\t-c, --cal-metafits=FILE    FILE is the metafits file pertaining to the calibration solution\n"
            "\t-C, --cal-location=PATH    PATH is the directory (RTS) or the file (OFFRINGA) containing the calibration solution\n"
            "\t-F, --flagged-tiles=FILE   FILE is a text file containing the TileNames of tiles to be flagged\n"
            "\t-N, --ncoarse_chans=NUM    NUM is the number of coarse channels to include\n"
            "\t-R, --ref-ant=TILENAME     Override the reference tile given in pq_phase_correction.txt for rotating the phases\n"
            "\t                           of the PP and QQ elements of the calibration solution. To turn off phase rotation\n"
            "\t                           altogether, set TILENAME=NONE.\n"
            "\t-U, --PQ-phase=PH,OFFS     Override the phase correction given in pq_phase_correction.txt. PH is given in rad/Hz\n"
            "\t                           and OFFS given in rad, such that, the QQ element of the calibration Jones matrix\n"
            "\t                           for frequency F (in Hz) is multiplied by\n"
            "\t                                exp(PH*F + OFFS)\n"
            "\t                           Setting PH = OFFS = 0 is equivalent to not performing any phase correction\n"
            "\t-X, --cross-terms          Retain the PQ and QP terms of the calibration solution [default: off]\n"
          );

    printf( "\nCALIBRATION OPTIONS (OFFRINGA)\n\n"
            "\t-O, --offringa             The calibration solution is in the Offringa format instead of\n"
            "\t                           the default RTS format. In this case, the argument to -C should\n" 
            "\t                           be the full path to the binary solution file.\n"
          );

    printf( "\nOTHER OPTIONS\n\n"
            "\t-h, --help                 Print this help and exit\n"
            "\t-V, --version              Print version number and exit\n\n"
          );
}



void make_plot_calibrate_parse_cmdline( int argc, char **argv,
        struct make_plot_calibrate_opts *opts )
{
    // Set defaults
    opts->ncoarse_chans        = -1;    // Number of coarse channels to include
    opts->metafits             = NULL;  // filename of the metafits file for the target observations
    opts->custom_flags         = NULL;  // filename of text file containing TileNames of tiles to be flagged
    opts->zero_string          = NULL;  // string to output if all matrix elements are identically zero
    opts->use_bandpass         = false; // use the Bandpass calibration solutions
    opts->metafits             = NULL;  // filename of the metafits file for the calibration observation
    opts->caldir               = NULL;  // The path to where the calibration solutions live
    opts->cal_type             = CAL_RTS;
    opts->ref_ant              = NULL;
    opts->keep_cross_terms     = false;
    opts->phase_offset         = 0.0;
    opts->phase_slope          = 0.0;
    opts->custom_pq_correction = false;

    if (argc > 1) {

        int c;
        while (1)
        {
            static struct option long_options[] = {
                {"bandpass",        no_argument,       0, 'B'},
                {"cal-location",    required_argument, 0, 'C'},
                {"cal-metafits",    required_argument, 0, 'c'},
                {"flagged-tiles",   required_argument, 0, 'F'},
                {"help",            no_argument,       0, 'h'},
                {"metafits",        required_argument, 0, 'm'},
                {"ncoarse_chans",   required_argument, 0, 'N'},
                {"offringa",        no_argument      , 0, 'O'},
                {"ref-ant",         required_argument, 0, 'R'},
                {"PQ-phase",        required_argument, 0, 'U'},
                {"version",         required_argument, 0, 'V'},
                {"cross-terms",     no_argument,       0, 'X'},
                {"zero-string",     required_argument, 0, 'z'},
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "Bc:C:F:hm:N:OR:U:VXz:",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c)
            {
                case 'B':
                    opts->use_bandpass = true;
                    break;
                case 'c':
                    opts->cal_metafits = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->cal_metafits, optarg );
                    break;
                case 'C':
                    opts->caldir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->caldir, optarg );
                    break;
                case 'F':
                    opts->custom_flags = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->custom_flags, optarg );
                    break;
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'N':
                    opts->ncoarse_chans = atoi(optarg);
                    break;
                case 'O':
                    opts->cal_type = CAL_OFFRINGA;
                    break;
                case 'R':
                    opts->ref_ant = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->ref_ant, optarg );
                    break;
                case 'U':
                    if (sscanf( optarg, "%lf,%lf", &(opts->phase_slope), &(opts->phase_offset) ) != 2)
                    {
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "cannot parse -U option (\"%s\") as \"FLOAT,FLOAT\"\n", optarg );
                        exit(EXIT_FAILURE);
                    }
                    opts->custom_pq_correction = true;
                    break;
                case 'V':
                    printf( "MWA Beamformer %s\n", VCSBEAM_VERSION);
                    exit(0);
                    break;
                case 'X':
                    opts->keep_cross_terms = true;
                    break;
                case 'z':
                    opts->zero_string = optarg;
                    break;
                default:
                    fprintf( stderr, "error: make_plot_calibrate_parse_cmdline: "
                                    "unrecognised option '%s'\n", optarg );
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


    // Check that all the required options were supplied
    assert( opts->metafits       != NULL );
    assert( opts->caldir         != NULL );
    assert( opts->cal_metafits   != NULL );
}
