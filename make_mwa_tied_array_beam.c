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
#include "jones.h"
#include "beam_psrfits.h"
#include "beam_vdif.h"
#include "metadata.h"
#include "form_beam.h"
#include "calibration.h"
#include "primary_beam.h"
#include "geometry.h"
#include "performance.h"
#include "ipfb.h"

#define MAX_COMMAND_LENGTH 1024

struct make_tied_array_beam_opts {
    char              *begin_str;     // Absolute or relative GPS time -- when to start beamforming
    unsigned long int  nseconds;      // How many seconds to process
    char              *pointings_file; // Name of file containing pointings (e.g. "hh:mm:ss dd:mm:ss")
    char              *datadir;       // The path to where the recombined data live
    char              *metafits;      // filename of the metafits file
    char              *coarse_chan_str;   // Absolute or relative coarse channel number
    int                ncoarse_chans; // How many coarse channels to process

    // Variables for MWA/VCS configuration
    char              *custom_flags;  // Use custom list for flagging antennas

    // Output options
    int                out_incoh;     // Default = PSRFITS (incoherent) output turned OFF
    int                out_coh;       // Default = PSRFITS (coherent)   output turned OFF
    int                out_vdif;      // Default = VDIF                 output turned OFF
    int                out_uvdif;     // Default = upsampled VDIF       output turned OFF
    int                out_bf;        // Default = beamform over all (non-flagged) antennas
    int                out_ant;       // The antenna number (0-127) to write out if out_bf = 0

    // Other options
    char              *synth_filter;  // Which synthesis filter to use
    int                out_summed;    // Default = output only Stokes I output turned OFF
    int                max_sec_per_file;    // Number of seconds per fits files
    float              gpu_mem;       // Default = -1.0. If -1.0 use all GPU mem
};

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/

void usage();
void make_tied_array_beam_parse_cmdline( int argc, char **argv, struct make_tied_array_beam_opts *opts, struct calibration *cal );
void parse_pointing_file( const char *filename, double **ras_hours, double **decs_degs, unsigned int *npointings );

/********
 * MAIN *
 ********/

int main(int argc, char **argv)
{
    // Parse command line arguments
    struct make_tied_array_beam_opts opts;
    struct calibration cal;           // Variables for calibration settings
    make_tied_array_beam_parse_cmdline( argc, argv, &opts, &cal );

    // Start a logger for output messages and time-keeping
    logger *log = create_logger( stdout );
    logger_add_stopwatch( log, "read" );
    logger_add_stopwatch( log, "delay" );
    logger_add_stopwatch( log, "calc" );
    logger_add_stopwatch( log, "write" );
    char log_message[MAX_COMMAND_LENGTH];

    // Create an mwalib metafits context and associated metadata
    logger_timed_message( log, "Creating metafits and voltage contexts via MWALIB" );

    char error_message[ERROR_MESSAGE_LEN];

    MetafitsContext  *obs_context  = NULL;
    MetafitsMetadata *obs_metadata = NULL;
    get_mwalib_metafits_metadata( opts.metafits, &obs_metadata, &obs_context );

    unsigned long int begin_gps = parse_begin_string( obs_metadata, opts.begin_str );
    uintptr_t begin_coarse_chan_idx = parse_coarse_chan_string( obs_metadata, opts.coarse_chan_str );

    VoltageMetadata  *vcs_metadata = NULL;
    VoltageContext   *vcs_context  = NULL;
    get_mwalib_voltage_metadata( &vcs_metadata, &vcs_context, &obs_metadata, obs_context,
            begin_gps, opts.nseconds, opts.datadir, begin_coarse_chan_idx, &opts.ncoarse_chans );

    MetafitsContext  *cal_context  = NULL;
    MetafitsMetadata *cal_metadata = NULL;
    get_mwalib_metafits_metadata( cal.metafits, &cal_metadata, &cal_context );

    uintptr_t ntimesteps = vcs_metadata->num_provided_timesteps;

    // Create some "shorthand" variables for code brevity
    uintptr_t nants          = obs_metadata->num_ants;
    uintptr_t nchans         = obs_metadata->num_volt_fine_chans_per_coarse;
    //int chan_width           = obs_metadata->volt_fine_chan_width_hz;
    uintptr_t npols          = obs_metadata->num_ant_pols;   // (X,Y)
    uintptr_t ninput         = obs_metadata->num_rf_inputs;
    uintptr_t outpol_coh     = 4;  // (I,Q,U,V)
    if ( opts.out_summed )
        outpol_coh           = 1;  // (I)
    const uintptr_t outpol_incoh = 1;  // ("I")
    unsigned int nsamples = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;

    int offset;
    unsigned int s;
    uintptr_t ch;
    uintptr_t pol;

    uintptr_t timestep;
    uintptr_t timestep_idx;
    uint64_t  gps_second;

    uintptr_t data_size = vcs_metadata->num_voltage_blocks_per_timestep * vcs_metadata->voltage_block_size_bytes;

    primary_beam pb;
    geometric_delays gdelays;

    // Start counting time from here (i.e. after parsing the command line)
    sprintf( log_message, "Reading pointings file %s", opts.pointings_file );
    logger_timed_message( log, log_message );

    // Parse input pointings
    double *ras_hours, *decs_degs;
    unsigned int npointing;
    unsigned int p;
    parse_pointing_file( opts.pointings_file, &ras_hours, &decs_degs, &npointing );

    // Allocate memory for various data products
    cuDoubleComplex  ****invJi                 = create_invJi( nants, nchans, npols );
    cuDoubleComplex  ****detected_beam         = create_detected_beam( npointing, 2*nsamples, nchans, npols );


    double invw = 1.0/nants;

    // Get pointing geometry information
    struct beam_geom beam_geom_vals[npointing];

    double mjd, sec_offset;
    mjd = obs_metadata->sched_start_mjd;
    for (p = 0; p < npointing; p++)
        calc_beam_geom( ras_hours[p], decs_degs[p], mjd, &beam_geom_vals[p] );

    // Create structures for the PFB filter coefficients
    int ntaps, fil_size = 0;
    double *coeffs = NULL;

    // If no synthesis filter was explicitly chosen, choose the LSQ12 filter
    if (!opts.synth_filter)
    {
        opts.synth_filter = (char *)malloc( 6 );
        strcpy( opts.synth_filter, "LSQ12" );
    }
    if (strcmp( opts.synth_filter, "LSQ12" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchans; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = LSQ12_FILTER_COEFFS; // I'll have to change the way these coefficients are stored
                                                   // in order to avoid this cumbersome loading procedure
        for (int i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else if (strcmp( opts.synth_filter, "MIRROR" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchans; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = MIRROR_FILTER_COEFFS;
        for (int i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else
    {
        fprintf( stderr, "error: unrecognised synthesis filter: %s\n",
                opts.synth_filter );
        exit(EXIT_FAILURE);
    }
    cuDoubleComplex *twiddles = roots_of_unity( nchans );

    // Adjust by the scaling that was introduced by the forward PFB,
    // along with any other scaling that I, Lord and Master of the inverse
    // PFB, feel is appropriate.
    double approx_filter_scale = 15.0/7.2; // 7.2 = 16384/117964.8
    for (int i = 0; i < fil_size; i++)
        coeffs[i] *= approx_filter_scale;

    /*********************
     * Memory Allocation *
     *********************/

    uint8_t *data;
    cudaMallocHost( (void **)&data, data_size );
    cudaCheckErrors( "cudaMallocHost(data) failed" );

    /* Allocate host and device memory for the use of the cu_form_beam function */
    // Declaring pointers to the structs so the memory can be alternated
    struct gpu_formbeam_arrays gf;

    struct gpu_ipfb_arrays gi;
    int nchunk;
    malloc_formbeam( &gf, nsamples, nants, nchans, npols, &nchunk, opts.gpu_mem,
                     outpol_coh, outpol_incoh, npointing, log );

    // Create a lists of rf_input indexes ordered by antenna number (needed for gpu kernels)
    create_antenna_lists( obs_metadata, gf.polX_idxs, gf.polY_idxs );

    // ... and upload them to the gpu, ready for use!
    cu_upload_pol_idx_lists( &gf );

    // Create output buffer arrays
    float *data_buffer_coh    = NULL;
    float *data_buffer_incoh  = NULL;
    float *data_buffer_vdif   = NULL;

    data_buffer_coh   = create_pinned_data_buffer_psrfits( npointing * nchans * outpol_coh * nsamples );
    data_buffer_incoh = create_pinned_data_buffer_psrfits( nchans * outpol_incoh * nsamples );
    data_buffer_vdif  = create_pinned_data_buffer_vdif( nsamples * nchans * npols * npointing * 2 * sizeof(float) );

    if (opts.out_vdif)
    {
        malloc_ipfb( &gi, ntaps, nsamples, nchans, npols, fil_size, npointing );
        cu_load_filter( coeffs, twiddles, &gi, nchans );
    }

    // Set up parallel streams
    cudaStream_t streams[npointing];
    for (p = 0; p < npointing; p++)
        cudaStreamCreate(&(streams[p])) ;

    // Create structures for holding header information
    struct psrfits  *pfs;
    struct psrfits   pf_incoh;
    pfs = (struct psrfits *)malloc(npointing * sizeof(struct psrfits));
    vdif_header     vhdr;
    struct vdifinfo *vf;
    vf = (struct vdifinfo *)malloc(npointing * sizeof(struct vdifinfo));

    cuDoubleComplex ***D = NULL; // See Eqs. (27) to (29) in Ord et al. (2019)

    // COARSE CHANNEL DEPENDENT CODE BEGINS HERE
    for (uintptr_t coarse_chan_idx = begin_coarse_chan_idx; coarse_chan_idx < begin_coarse_chan_idx + opts.ncoarse_chans; coarse_chan_idx++)
    {
        /****************************
         * GET CALIBRATION SOLUTION *
         ****************************/

        if (cal.cal_type == CAL_RTS)
        {
            D = get_rts_solution( cal_metadata, obs_metadata, cal.caldir, coarse_chan_idx );
            if (cal.apply_xy_correction)
            {
                xy_phase_correction( obs_metadata->obs_id, &cal.phase_slope, &cal.phase_offset );

                // Print a suitable message
                if (cal.phase_slope == 0.0 && cal.phase_offset == 0.0)
                    logger_timed_message( log, "No XY phase correction information for this obsid" );
                else
                {
                    sprintf( log_message, "Applying XY phase correction %.2e*freq%+.2e",
                            cal.phase_slope, cal.phase_offset );
                    logger_timed_message( log, log_message );
                }
            }
            else
            {
                logger_timed_message( log, "Not applying XY phase correction" );
            }
        }
        else if (cal.cal_type == CAL_OFFRINGA)
        {
            fprintf( stderr, "error: Offringa-style calibration solutions not currently supported\n" );
            exit(EXIT_FAILURE);
            /*
            // Find the ordering of antennas in Offringa solutions from metafits file
            read_offringa_gains_file( D, nants, cal.offr_chan_num, cal.filename );
            */
        }

        // ------------------
        // Prepare primary beam and geometric delay arrays
        // ------------------
        create_primary_beam( &pb, obs_metadata, coarse_chan_idx, npointing );
        create_geometric_delays( &gdelays, obs_metadata, vcs_metadata, coarse_chan_idx, npointing );

        // ------------------

        sprintf( log_message, "Preparing headers for output (receiver channel %lu)",
                obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
        logger_message( log, log_message );

        // Populate the relevant header structs
        for (p = 0; p < npointing; p++)
            populate_psrfits_header( &pfs[p], obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
                    outpol_coh, &beam_geom_vals[p], NULL, true );

        populate_psrfits_header( &pf_incoh, obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
                outpol_incoh, &beam_geom_vals[0], NULL, false );

        populate_vdif_header( vf, &vhdr, obs_metadata, vcs_metadata, coarse_chan_idx,
                beam_geom_vals, npointing );

        // Begin the main loop: go through data one second at a time

        logger_message( log, "\n*****BEGIN BEAMFORMING*****" );

        for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
        {
            logger_message( log, "" ); // Print a blank line

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
                fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s", error_message );
                exit(EXIT_FAILURE);
            }
            logger_stop_stopwatch( log, "read" );

            // Get the next second's worth of phases / jones matrices, if needed
            sprintf( log_message, "[%lu/%lu] Calculating Jones matrices",
                    timestep_idx+1, ntimesteps );
            logger_timed_message( log, log_message );

            logger_start_stopwatch( log, "delay" );
            sec_offset = (double)(timestep_idx + begin_gps - obs_metadata->obs_id);
            mjd = obs_metadata->sched_start_mjd + (sec_offset + 0.5)/86400.0;
            for (p = 0; p < npointing; p++)
                calc_beam_geom( ras_hours[p], decs_degs[p], mjd, &beam_geom_vals[p] );

            // Calculate the primary beam
            calc_primary_beam( &pb, beam_geom_vals );

            // Calculate the geometric delays
            calc_all_geometric_delays( &gdelays, beam_geom_vals );
            push_geometric_delays_to_device( &gdelays );

            get_jones(
                    npointing,              // number of pointings
                    obs_metadata,
                    coarse_chan_idx,
                    &cal,                   // struct holding info about calibration
                    D,                      // Calibration Jones matrices
                    pb.B,                   // Primary beam jones matrices
                    invJi );                // invJi array           (answer will be output here)

            logger_stop_stopwatch( log, "delay" );

            // Form the beams
            sprintf( log_message, "[%lu/%lu] Calculating beam", timestep_idx+1, ntimesteps);
            logger_timed_message( log, log_message );

            logger_start_stopwatch( log, "calc" );
            if (!opts.out_bf) // don't beamform, but only procoess one ant/pol combination
            {
                // Populate the detected_beam, data_buffer_coh, and data_buffer_incoh arrays
                // detected_beam = [2*nsamples][nchans][npols] = [20000][128][2] (cuDoubleComplex)
                if (timestep_idx % 2 == 0)
                    offset = 0;
                else
                    offset = nsamples;

                for (p   = 0; p   < npointing;   p++  )
                    for (s   = 0; s   < nsamples; s++  )
                        for (ch  = 0; ch  < nchans;       ch++ )
                            for (pol = 0; pol < npols;        pol++)
                            {
                                int i;
                                if (pol == 0)
                                    i = gf.polX_idxs[opts.out_ant];
                                else
                                    i = gf.polY_idxs[opts.out_ant];
                                detected_beam[p][s+offset][ch][pol] = UCMPLX4_TO_CMPLX_FLT(data[v_IDX(s,ch,i,nchans,ninput)]);
                                detected_beam[p][s+offset][ch][pol] = cuCmul(detected_beam[p][s+offset][ch][pol], make_cuDoubleComplex(128.0, 0.0));
                            }
            }
            else // beamform (the default mode)
            {
                cu_form_beam( data, nsamples, gdelays.d_phi, invJi, timestep_idx,
                        npointing, nants, nchans, npols, outpol_coh, invw, &gf,
                        detected_beam, data_buffer_coh, data_buffer_incoh,
                        streams, opts.out_incoh, nchunk );
            }

            // Invert the PFB, if requested
            if (opts.out_vdif)
            {
                sprintf( log_message, "[%lu/%lu] Inverting the PFB (full)",
                        timestep_idx+1, ntimesteps);
                logger_timed_message( log, log_message );
                cu_invert_pfb_ord( detected_beam, timestep_idx, npointing,
                        nsamples, nchans, npols, vf->sizeof_buffer,
                        &gi, data_buffer_vdif );
            }
            logger_stop_stopwatch( log, "calc" );


            // Write out for each pointing
            for (p = 0; p < npointing; p++)
            {
                sprintf( log_message, "[%lu/%lu] [%d/%d] Writing data to file(s)",
                        timestep_idx+1, ntimesteps, p+1, npointing );
                logger_timed_message( log, log_message );

                logger_start_stopwatch( log, "write" );

                if (opts.out_coh)
                    psrfits_write_second( &pfs[p], data_buffer_coh, nchans,
                            outpol_coh, p );
                if (opts.out_incoh && p == 0)
                    psrfits_write_second( &pf_incoh, data_buffer_incoh,
                            nchans, outpol_incoh, p );
                if (opts.out_vdif)
                    vdif_write_second( &vf[p], &vhdr,
                            data_buffer_vdif + p * vf->sizeof_buffer );
                logger_stop_stopwatch( log, "write" );
            }
        }

        logger_message( log, "\n*****END BEAMFORMING*****\n" );

        // Clean up channel-dependent memory
        if (opts.out_incoh)
        {
            free_psrfits( &pf_incoh );
        }
        for (p = 0; p < npointing; p++)
        {
            if (opts.out_coh)
            {
                free_psrfits( &pfs[p] );
            }
            if (opts.out_vdif)
            {
                free( vf[p].b_scales  );
                free( vf[p].b_offsets );
            }
        }

        free_rts( D, cal_metadata );
    }

    logger_report_all_stats( log );
    logger_message( log, "" );

    // Free up memory
    logger_timed_message( log, "Starting clean-up" );

    destroy_invJi( invJi, nants, nchans, npols );
    destroy_detected_beam( detected_beam, npointing, 2*nsamples, nchans );

    free( twiddles );
    free( coeffs );

    cudaFreeHost( data_buffer_coh   );
    cudaCheckErrors( "cudaFreeHost(data_buffer_coh) failed" );
    cudaFreeHost( data_buffer_incoh );
    cudaCheckErrors( "cudaFreeHost(data_buffer_incoh) failed" );
    cudaFreeHost( data_buffer_vdif  );
    cudaCheckErrors( "cudaFreeHost(data_buffer_vdif) failed" );
    cudaFreeHost( data );
    cudaCheckErrors( "cudaFreeHost(data) failed" );

    free( opts.pointings_file  );
    free( opts.datadir         );
    free( opts.begin_str       );
    free( opts.coarse_chan_str );
    free( cal.caldir           );
    free( opts.metafits        );
    free( opts.synth_filter    );

    free_formbeam( &gf );
    if (opts.out_vdif)
    {
        free_ipfb( &gi );
    }

    // Clean up memory associated with the Jones matrices
    free_primary_beam( &pb );
    free_geometric_delays( &gdelays );

    // Clean up memory associated with mwalib
    mwalib_metafits_metadata_free( obs_metadata );
    mwalib_voltage_metadata_free( vcs_metadata );
    mwalib_voltage_context_free( vcs_context );
    mwalib_metafits_metadata_free( cal_metadata );
    mwalib_metafits_context_free( cal_context );


    return EXIT_SUCCESS;
}


void usage() {
    printf( "\nusage: make_mwa_tied_array_beam [OPTIONS]\n");

    printf( "\nREQUIRED OPTIONS\n\n"
            "\t-P, --pointings=FILE       FILE containing RA and Decs of multiple pointings\n"
            "\t                           in the format hh:mm:ss.s dd:mm:ss.s ...\n"
            "\t-m, --metafits=FILE        FILE is the metafits file for the target observation\n"
          );

    printf( "\nOPTIONAL OPTIONS\n\n"
            "\t-b, --begin=GPSTIME        Begin time of observation, in GPS seconds\n"
            "\t                           If GPSTIME starts with a '+' or a '-', then the time\n"
            "\t                           is taken relative to the start or end of the observation\n"
            "\t                           respectively. [default: \"+0\"]\n"
            "\t-d, --data-location=PATH   PATH is the directory containing the recombined data\n"
            "\t                           [default: current directory]\n"
            "\t-f, --coarse-chan=CHAN     Coarse channel number\n"
            "\t                           If CHAN starts with a '+' or a '-', then the channel is taken\n"
            "\t                           relative to the first or last channel in the observation\n"
            "\t                           respectively. Otherwise, it is treated as a receiver channel number\n"
            "\t                           (0-255) [default: \"+0\"]\n"
            "\t-F, --nchans=VAL           Process VAL coarse channels\n"
            "\t-i, --incoh                Turn on incoherent PSRFITS beam output.                             [default: OFF]\n"
            "\t-p, --psrfits              Turn on coherent PSRFITS output (will be turned on if none of\n"
            "\t                           -i, -p, -u, -v are chosen).                                         [default: OFF]\n"
            "\t-v, --vdif                 Turn on VDIF output with upsampling                                 [default: OFF]\n"
            "\t-s, --summed               Turn on summed polarisations of the coherent output (only Stokes I) [default: OFF]\n"
            "\t-t, --max_t                Maximum number of seconds per output fits file. [default: 200]\n"
            "\t-A, --antpol=ant           Do not beamform. Instead, only operate on the specified ant\n"
            "\t                           stream (0-127)\n" 
            "\t-S, --synth_filter=filter  Apply the named filter during high-time resolution synthesis.\n"
            "\t                           filter can be MIRROR or LSQ12.\n"
            "\t                           [default: LSQ12]\n"
            "\t-T, --nseconds=VAL         Process VAL seconds of data [default: as many as possible]\n"
          );

    printf( "\nCALIBRATION OPTIONS (RTS)\n\n"
            "\t-c, --cal-metafits=FILE    FILE is the metafits file pertaining to the calibration solution\n"
            "\t-C, --cal-location=PATH    PATH is the directory (RTS) or the file (OFFRINGA) containing the calibration solution\n"
            "\t-R, --ref-ant=ANT          Rotate the phases of the XX and YY elements of the calibration\n"
            "\t                           Jones matrices so that the phases of tile ANT align. If ANT is\n"
            "\t                           outside the range 0-127, no phase rotation is done\n"
            "\t                           [default: 0]\n"
            "\t-X, --cross-terms          Retain the XY and YX terms of the calibration solution\n"
            "\t                           [default: off]\n"
            "\t-U, --no-XY-phase          Do not apply the XY phase correction to the calibration solution\n"
          );

    printf( "\nCALIBRATION OPTIONS (OFFRINGA) -- NOT YET SUPPORTED\n\n"
            "\t-O, --offringa             The calibration solution is in the Offringa format instead of\n"
            "\t                           the default RTS format. In this case, the argument to -C should\n" 
            "\t                           be the full path to the binary solution file.\n"
          );

    printf( "\nOTHER OPTIONS\n\n"
            "\t-h, --help                 Print this help and exit\n"
            "\t-g, --gpu-mem=N            The maximum amount of GPU memory you want make_beam to use in GB [default: -1]\n"
            "\t-V, --version              Print version number and exit\n\n"
          );
}



void make_tied_array_beam_parse_cmdline(
        int argc, char **argv, struct make_tied_array_beam_opts *opts, struct calibration *cal )
{
    // Set defaults
    opts->begin_str          = NULL; // Absolute or relative GPS time -- when to start beamforming
    opts->nseconds           = -1;   // How many seconds to process (-1 = as many as possible)
    opts->pointings_file     = NULL; // File containing list of pointings "hh:mm:ss dd:mm:ss ..."
    opts->datadir            = NULL; // The path to where the recombined data live
    opts->metafits           = NULL; // filename of the metafits file for the target observation
    opts->coarse_chan_str    = NULL; // Absolute or relative coarse channel
    opts->ncoarse_chans      = -1;   // How many coarse channels to process
    opts->out_incoh          = 0;    // Default = PSRFITS (incoherent) output turned OFF
    opts->out_coh            = 0;    // Default = PSRFITS (coherent)   output turned OFF
    opts->out_vdif           = 0;    // Default = VDIF                 output turned OFF
    opts->out_bf             = 1;    // Default = beamform all (non-flagged) antennas
    opts->out_ant            = 0;    // The antenna number (0-127) to write out if out_bf = 0
    opts->synth_filter       = NULL;
    opts->out_summed         = 0;    // Default = output only Stokes I output turned OFF
    opts->max_sec_per_file   = 200;  // Number of seconds per fits files
    opts->gpu_mem            = -1.0;

    cal->metafits            = NULL; // filename of the metafits file for the calibration observation
    cal->caldir              = NULL; // The path to where the calibration solutions live
    cal->cal_type            = CAL_RTS;
    cal->ref_ant             = 0;
    cal->cross_terms         = 0;
    cal->phase_offset        = 0.0;
    cal->phase_slope         = 0.0;
    cal->apply_xy_correction = true;

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"incoh",           no_argument,       0, 'i'},
                {"psrfits",         no_argument,       0, 'p'},
                {"vdif",            no_argument,       0, 'v'},
                {"summed",          no_argument,       0, 's'},
                {"max_t",           required_argument, 0, 't'},
                {"synth_filter",    required_argument, 0, 'S'},
                {"nseconds",        required_argument, 0, 'T'},
                {"antpol",          required_argument, 0, 'A'},
                {"pointings",       required_argument, 0, 'P'},
                {"data-location",   required_argument, 0, 'd'},
                {"cal-location",    required_argument, 0, 'C'},
                {"metafits",        required_argument, 0, 'm'},
                {"cal-metafits",    required_argument, 0, 'c'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"nchans",          required_argument, 0, 'F'},
                {"ref-ant",         required_argument, 0, 'R'},
                {"cross-terms",     no_argument,       0, 'X'},
                {"no-XY-phase",     required_argument, 0, 'U'},
                {"offringa",        no_argument      , 0, 'O'},
                {"gpu-mem",         required_argument, 0, 'g'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "A:b:c:d:C:f:F:g:him:OpP:R:sS:t:T:UvVX",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'A':
                    opts->out_bf = 0; // Turn off normal beamforming
                    opts->out_ant = atoi(optarg); // 0-127
                    break;
                case 'b':
                    opts->begin_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->begin_str, optarg );
                    break;
                case 'c':
                    cal->metafits = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( cal->metafits, optarg );
                    break;
                case 'C':
                    cal->caldir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( cal->caldir, optarg );
                    break;
                case 'd':
                    opts->datadir = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->datadir, optarg );
                    break;
                case 'f':
                    opts->coarse_chan_str = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->coarse_chan_str, optarg );
                    break;
                case 'F':
                    opts->ncoarse_chans = atoi(optarg);
                    break;
                case 'g':
                    opts->gpu_mem = atof(optarg);
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'i':
                    opts->out_incoh = 1;
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'O':
                    cal->cal_type = CAL_OFFRINGA;
                    break;
                case 'p':
                    opts->out_coh = 1;
                    break;
                case 'P':
                    opts->pointings_file = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->pointings_file, optarg );
                    break;
                case 'R':
                    cal->ref_ant = atoi(optarg);
                    break;
                case 'S':
                    opts->synth_filter = (char *)malloc( strlen(optarg) + 1 );
                    strcpy( opts->synth_filter, optarg );
                    break;
                case 's':
                    opts->out_summed = 1;
                    break;
                case 't':
                    opts->max_sec_per_file = atoi(optarg);
                    break;
                case 'T':
                    opts->nseconds = atol(optarg);
                    if (opts->nseconds <= 0)
                    {
                        fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
                                "-%c argument must be >= 1\n", c );
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'U':
                    cal->apply_xy_correction = false;
                    break;
                case 'v':
                    opts->out_vdif = 1;
                    break;
                case 'V':
                    printf( "MWA Beamformer %s\n", VERSION_BEAMFORMER);
                    exit(0);
                    break;
                case 'X':
                    cal->cross_terms = 1;
                    break;
                default:
                    fprintf( stderr, "error: make_tied_array_beam_parse_cmdline: "
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
    assert( opts->pointings_file != NULL );
    assert( cal->caldir          != NULL );
    assert( opts->metafits       != NULL );
    assert( cal->metafits        != NULL );

    if (opts->datadir == NULL)
    {
        opts->datadir = (char *)malloc( 2 );
        strcpy( opts->datadir, "." );
    }

    if (opts->begin_str == NULL)
    {
        opts->begin_str = (char *)malloc( 3 );
        strcpy( opts->begin_str, "+0" );
    }

    if (opts->coarse_chan_str == NULL)
    {
        opts->coarse_chan_str = (char *)malloc( 3 );
        strcpy( opts->coarse_chan_str, "+0" );
    }

    // If neither -i, -p, nor -v were chosen, set -p by default
    if ( !opts->out_incoh && !opts->out_coh && !opts->out_vdif )
    {
        opts->out_coh = 1;
    }

}



cuDoubleComplex ****create_invJi( int nants, int nchans, int npols )
// Allocate memory for (inverse) Jones matrices
{
    int ant, pol, ch; // Loop variables
    cuDoubleComplex ****invJi;
    invJi = (cuDoubleComplex ****)malloc( nants * sizeof(cuDoubleComplex ***) );

    for (ant = 0; ant < nants; ant++)
    {
        invJi[ant] =(cuDoubleComplex ***)malloc( nchans * sizeof(cuDoubleComplex **) );

        for (ch = 0; ch < nchans; ch++)
        {
            invJi[ant][ch] = (cuDoubleComplex **)malloc( npols * sizeof(cuDoubleComplex *) );

            for (pol = 0; pol < npols; pol++)
                invJi[ant][ch][pol] = (cuDoubleComplex *)malloc( npols * sizeof(cuDoubleComplex) );
        }
    }
    return invJi;
}


void destroy_invJi( cuDoubleComplex ****array, int nants, int nchans, int npols )
{
    int ant, ch, pol;
    for (ant = 0; ant < nants; ant++)
    {
        for (ch = 0; ch < nchans; ch++)
        {
            for (pol = 0; pol < npols; pol++)
                free( array[ant][ch][pol] );

            free( array[ant][ch] );
        }
        free( array[ant] );
    }
    free( array );
}


void parse_pointing_file( const char *filename, double **ras_hours, double **decs_degs, unsigned int *npointings )
/* Parse the given file in FILENAME and create arrays of RAs and DECs.
 * This function allocates memory for ras_hours and decs_degs arrays.
 * Caller can destroy with free().
 */
{
    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error: cannot open pointings file %s\n", filename );
        exit(EXIT_FAILURE);
    }

    // Do one pass through the file to count "words"
    // The RAs and Decs are expected to be whitespace-delimited
    int nwords = 0;
    char word[64];
    while (fscanf( f, "%s", word ) != EOF)
        nwords++;

    // Check that we have an even number of words (they should be in RA/Dec pairs)
    if (nwords % 2 != 0)
    {
        fprintf( stderr, "error: cannot parse pointings file %s\n", filename );
        exit(EXIT_FAILURE);
    }
    *npointings = nwords/2;

    // Allocate memory
    *ras_hours = (double *)malloc( *npointings * sizeof(double) );
    *decs_degs = (double *)malloc( *npointings * sizeof(double) );

    // Rewind to beginning of file, read the words in again, and parse them
    rewind( f );
    char ra_str[64], dec_str[64];
    unsigned int p;
    for (p = 0; p < *npointings; p++)
    {
        // Read in the next Ra/Dec pair
        fscanf( f, "%s %s", ra_str, dec_str );

        // Parse them and make them decimal
        (*ras_hours)[p] = parse_ra( ra_str );
        (*decs_degs)[p] = parse_dec( dec_str );
    }

    // Close the file
    fclose( f );
}
