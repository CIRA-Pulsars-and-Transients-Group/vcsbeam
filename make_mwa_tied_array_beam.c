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
#include "beam_vdif.h"
#include "metadata.h"
#include "vdifio.h"
#include "filter.h"
#include "psrfits.h"
#include "form_beam.h"
#include "calibration.h"
#include "primary_beam.h"
#include "geometric_delay.h"

#include <cuda_runtime.h>
#include "ipfb.h"

double now()
{
  struct timespec t;
  clock_gettime(CLOCK_REALTIME,&t);
  return (double)t.tv_sec + (double)t.tv_nsec/1000000000L;
}

#define NOW now()
#define MAX_COMMAND_LENGTH 1024

struct make_beam_opts {
    // Variables for required options
    unsigned long int  begin;         // GPS time -- when to start beamforming
    unsigned long int  end;           // GPS time -- when to stop beamforming
    char              *pointings_file; // Name of file containing pointings (e.g. "hh:mm:ss dd:mm:ss")
    char              *datadir;       // The path to where the recombined data live
    char              *metafits;      // filename of the metafits file
    uintptr_t          rec_channel;   // 0 - 255 receiver 1.28MHz channel
    long int           frequency;     // = rec_channel expressed in Hz

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
    float              gpu_mem  ;     // Default = -1.0. If -1.0 use all GPU mem
};

void usage();
void make_beam_parse_cmdline( int argc, char **argv, struct make_beam_opts *opts, struct calibration *cal );

int main(int argc, char **argv)
{
    // A place to hold the beamformer settings
    struct make_beam_opts opts;
    struct calibration cal;           // Variables for calibration settings

    int i; // Generic counter

    /* Set default beamformer settings */

    // Variables for required options
    opts.begin       = 0;    // GPS time -- when to start beamforming
    opts.end         = 0;    // GPS time -- when to stop beamforming
    opts.pointings_file = NULL; // File containing list of pointings "hh:mm:ss dd:mm:ss ..."
    opts.datadir     = NULL; // The path to where the recombined data live
    opts.metafits    = NULL; // filename of the metafits file for the target observation
    opts.rec_channel = -1;   // 0 - 255 receiver 1.28MHz channel

    // Output options
    opts.out_incoh     = 0;  // Default = PSRFITS (incoherent) output turned OFF
    opts.out_coh       = 0;  // Default = PSRFITS (coherent)   output turned OFF
    opts.out_vdif      = 0;  // Default = VDIF                 output turned OFF
    opts.out_bf        = 1;  // Default = beamform all (non-flagged) antennas
    opts.out_ant       = 0;  // The antenna number (0-127) to write out if out_bf = 0
    opts.synth_filter  = NULL;
    opts.out_summed    = 0;  // Default = output only Stokes I output turned OFF
    opts.max_sec_per_file = 200; // Number of seconds per fits files

    // Variables for calibration settings
    cal.metafits          = NULL; // filename of the metafits file for the calibration observation
    cal.caldir            = NULL; // The path to where the calibration solutions live
    cal.cal_type          = CAL_NONE;
    cal.ref_ant           = 0;
    cal.cross_terms       = 0;
    cal.phase_offset      = 0.0;
    cal.phase_slope       = 0.0;
    cal.apply_xy_correction = true;

    // GPU options
    opts.gpu_mem          = -1.0;

    // Parse command line arguments
    make_beam_parse_cmdline( argc, argv, &opts, &cal );


    double begintime = NOW;

    char error_message[ERROR_MESSAGE_LEN];

    FILE *log = stdout;

    // Create an mwalib metafits context and associated metadata
    fprintf( log, "[%f]  Creating metafits and voltage contexts via MWALIB\n", NOW-begintime );

    MetafitsMetadata *obs_metadata = NULL;
    MetafitsMetadata *cal_metadata = NULL;
    VoltageMetadata  *vcs_metadata = NULL;
    VoltageContext   *vcs_context  = NULL;

    get_mwalib_metadata( &obs_metadata, &vcs_metadata, &vcs_context, &cal_metadata,
            opts.metafits, cal.metafits, opts.begin, opts.end, opts.datadir, opts.rec_channel );

    uintptr_t ntimesteps = vcs_metadata->num_common_timesteps;

    // Create some "shorthand" variables for code brevity
    uintptr_t nant           = obs_metadata->num_ants;
    uintptr_t nchan          = obs_metadata->num_volt_fine_chans_per_coarse;
    //int chan_width           = obs_metadata->volt_fine_chan_width_hz;
    uintptr_t npol           = obs_metadata->num_ant_pols;   // (X,Y)
    uintptr_t ninput         = obs_metadata->num_rf_inputs;
    uintptr_t outpol_coh     = 4;  // (I,Q,U,V)
    if ( opts.out_summed )
        outpol_coh           = 1;  // (I)
    const uintptr_t outpol_incoh   = 1;  // ("I")
    unsigned int sample_rate = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;

    int offset;
    unsigned int s;
    uintptr_t ch;
    uintptr_t pol;


    // Start counting time from here (i.e. after parsing the command line)
    fprintf( log, "[%f]  Reading pointings file %s\n", NOW-begintime, opts.pointings_file );

    // Parse input pointings
    double *ras_hours, *decs_degs;
    unsigned int npointing;
    unsigned int p;
    parse_pointing_file( opts.pointings_file, &ras_hours, &decs_degs, &npointing );

    // Allocate memory for various data products
    cuDoubleComplex  ****invJi                 = create_invJi( nant, nchan, npol );
    cuDoubleComplex  ****detected_beam         = create_detected_beam( npointing, 2*sample_rate, nchan, npol );


    // If using bandpass calibration solutions, calculate number of expected bandpass channels
    double invw = 1.0/nant;

    // ------------------------
    // GET CALIBRATION SOLUTION
    // ------------------------
    cuDoubleComplex ***D = NULL; // See Eqs. (27) to (29) in Ord et al. (2019)
    if (cal.cal_type == CAL_RTS)
    {
        D = get_rts_solution( cal_metadata, obs_metadata, cal.caldir, opts.rec_channel );
        if (cal.apply_xy_correction)
        {
            xy_phase_correction( obs_metadata->obs_id, &cal.phase_slope, &cal.phase_offset );

            // Print a suitable message
            if (cal.phase_slope == 0.0 && cal.phase_offset == 0.0)
                fprintf( log, "[%f]  No XY phase correction information for this obsid\n",
                        NOW-begintime );
            else
                fprintf( log, "[%f]  Applying XY phase correction %.2e*freq%+.2e\n",
                        NOW-begintime, cal.phase_slope, cal.phase_offset );
        }
        else
        {
            fprintf( log, "[%f]  Not applying XY phase correction\n", NOW-begintime );
        }
    }
    else if (cal.cal_type == CAL_OFFRINGA)
    {
        fprintf( stderr, "error: Offringa-style calibration solutions not currently supported\n" );
        exit(EXIT_FAILURE);
        /*
        // Find the ordering of antennas in Offringa solutions from metafits file
        read_offringa_gains_file( D, nant, cal.offr_chan_num, cal.filename );
        */
    }

    int coarse_chan_idx = 0; /* Value is fixed for now (i.e. each call of make_beam only
                                ever processes one coarse chan. However, in the future,
                                this should be flexible, with mpi or threads managing
                                different coarse channels. */
    int coarse_chan = vcs_metadata->provided_coarse_chan_indices[coarse_chan_idx];

    // ------------------
    // Allocate memory for various arrays
    // ------------------
    primary_beam pb;
    create_primary_beam( &pb, obs_metadata, coarse_chan, npointing );

    geometric_delays gdelays;
    create_geometric_delays( &gdelays, obs_metadata, vcs_metadata, coarse_chan, npointing );

    // ------------------

    fprintf( log, "[%f]  Setting up output header information\n", NOW-begintime);
    struct beam_geom beam_geom_vals[npointing];

    double mjd, sec_offset;
    mjd = obs_metadata->sched_start_mjd;
    calc_beam_geom( ras_hours, decs_degs, npointing, mjd, beam_geom_vals );

    // Create structures for holding header information
    struct psrfits  *pf;
    struct psrfits  *pf_incoh;
    pf = (struct psrfits *)malloc(npointing * sizeof(struct psrfits));
    pf_incoh = (struct psrfits *)malloc(1 * sizeof(struct psrfits));
    vdif_header     vhdr;
    struct vdifinfo *vf;
    vf = (struct vdifinfo *)malloc(npointing * sizeof(struct vdifinfo));


    // Create structures for the PFB filter coefficients
    int ntaps, fil_size = 0;
    double *coeffs = NULL;

    // If no synthesis filter was explicitly chosen, choose the LSQ12 filter
    if (!opts.synth_filter)  opts.synth_filter = strdup("LSQ12");
    if (strcmp( opts.synth_filter, "LSQ12" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchan; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = LSQ12_FILTER_COEFFS; // I'll have to change the way these coefficients are stored
                                                   // in order to avoid this cumbersome loading procedure
        for (i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else if (strcmp( opts.synth_filter, "MIRROR" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchan; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = MIRROR_FILTER_COEFFS;
        for (i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else
    {
        fprintf( stderr, "error: unrecognised synthesis filter: %s\n",
                opts.synth_filter );
        exit(EXIT_FAILURE);
    }
    cuDoubleComplex *twiddles = roots_of_unity( nchan );

    // Adjust by the scaling that was introduced by the forward PFB,
    // along with any other scaling that I, Lord and Master of the inverse
    // PFB, feels is appropriate.
    double approx_filter_scale = 15.0/7.2; // 7.2 = 16384/117964.8
    for (i = 0; i < fil_size; i++)
        coeffs[i] *= approx_filter_scale;

    // Populate the relevant header structs
    populate_psrfits_header( pf, obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
            outpol_coh, beam_geom_vals, npointing, true );

    populate_psrfits_header( pf_incoh, obs_metadata, vcs_metadata, coarse_chan_idx, opts.max_sec_per_file,
            outpol_incoh, beam_geom_vals, 1, false );

    populate_vdif_header( vf, &vhdr, obs_metadata, vcs_metadata, coarse_chan_idx,
            beam_geom_vals, npointing );

    // To run asynchronously we require two memory allocations for each data
    // set so multiple parts of the memory can be worked on at once.
    // We control this by changing the pointer to alternate between
    // the two memory allocations

    // Create array for holding the raw data
    long int bytes_per_timestep = vcs_metadata->num_voltage_blocks_per_timestep * vcs_metadata->voltage_block_size_bytes;

    uint8_t *data;
    cudaMallocHost( (void **)&data, bytes_per_timestep * sizeof(uint8_t) );
    cudaCheckErrors( "cudaMallocHost(data) failed" );

    /* Allocate host and device memory for the use of the cu_form_beam function */
    // Declaring pointers to the structs so the memory can be alternated
    struct gpu_formbeam_arrays gf;

    struct gpu_ipfb_arrays gi;
    int nchunk;
    malloc_formbeam( &gf, sample_rate, nant, nchan, npol, &nchunk, opts.gpu_mem,
                     outpol_coh, outpol_incoh, npointing, NOW-begintime );

    // Create a lists of rf_input indexes ordered by antenna number (needed for gpu kernels)
    create_antenna_lists( obs_metadata, gf.polX_idxs, gf.polY_idxs );

    // ... and upload them to the gpu, ready for use!
    cu_upload_pol_idx_lists( &gf );

    // Create output buffer arrays
    float *data_buffer_coh    = NULL;
    float *data_buffer_incoh  = NULL;
    float *data_buffer_vdif   = NULL;

    data_buffer_coh   = create_pinned_data_buffer_psrfits( npointing * nchan *
                                                           outpol_coh * pf[0].hdr.nsblk );
    data_buffer_incoh = create_pinned_data_buffer_psrfits( nchan * outpol_incoh *
                                                           pf_incoh[0].hdr.nsblk );
    data_buffer_vdif  = create_pinned_data_buffer_vdif( vf->sizeof_buffer *
                                                        npointing );

    if (opts.out_vdif)
    {
        malloc_ipfb( &gi, ntaps, sample_rate, nchan, npol, fil_size, npointing );
        cu_load_filter( coeffs, twiddles, &gi, nchan );
    }

    // Set up parallel streams
    cudaStream_t streams[npointing];
    for (p = 0; p < npointing; p++)
        cudaStreamCreate(&(streams[p])) ;

    fprintf( log, "\n[%f]  *****BEGINNING BEAMFORMING*****\n", NOW-begintime );

    // Set up timing for each section
    long read_time[ntimesteps], delay_time[ntimesteps], calc_time[ntimesteps], write_time[ntimesteps][npointing];
    uintptr_t timestep_idx;
    long timestep;
    uint64_t gps_second;

    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        timestep = vcs_metadata->common_timestep_indices[timestep_idx];
        coarse_chan = vcs_metadata->provided_coarse_chan_indices[coarse_chan_idx];
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
                    bytes_per_timestep,
                    error_message,
                    ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
        {
            fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s\n", error_message );
            exit(EXIT_FAILURE);
        }
        read_time[timestep_idx] = clock() - start;

        // Get the next second's worth of phases / jones matrices, if needed
        start = clock();
        fprintf( log, "[%f] [%lu/%lu] Calculating Jones matrices\n", NOW-begintime,
                                timestep_idx+1, ntimesteps );

        sec_offset = (double)(timestep_idx + opts.begin - obs_metadata->obs_id);
        mjd = obs_metadata->sched_start_mjd + (sec_offset + 0.5)/86400.0;
        calc_beam_geom( ras_hours, decs_degs, npointing, mjd, beam_geom_vals );

        // Calculate the primary beam
        calc_primary_beam( &pb, beam_geom_vals );

        // Calculate the geometric delays
        calc_geometric_delays( &gdelays, beam_geom_vals );
        push_geometric_delays_to_device( &gdelays );

        get_jones(
                npointing,              // number of pointings
                vcs_metadata,
                obs_metadata,
                coarse_chan_idx,
                &cal,                   // struct holding info about calibration
                D,                      // Calibration Jones matrices
                pb.B,                   // Primary beam jones matrices
                invJi );                // invJi array           (answer will be output here)

        delay_time[timestep_idx] = clock() - start;

        fprintf( log, "[%f] [%lu/%lu] Calculating beam\n", NOW-begintime,
                                timestep_idx+1, ntimesteps);
        start = clock();

        if (!opts.out_bf) // don't beamform, but only procoess one ant/pol combination
        {
            // Populate the detected_beam, data_buffer_coh, and data_buffer_incoh arrays
            // detected_beam = [2*sample_rate][nchan][npol] = [20000][128][2] (cuDoubleComplex)
            if (timestep_idx % 2 == 0)
                offset = 0;
            else
                offset = sample_rate;

            for (p   = 0; p   < npointing;   p++  )
            for (s   = 0; s   < sample_rate; s++  )
            for (ch  = 0; ch  < nchan;       ch++ )
            for (pol = 0; pol < npol;        pol++)
            {
                if (pol == 0)
                    i = gf.polX_idxs[opts.out_ant];
                else
                    i = gf.polY_idxs[opts.out_ant];
                detected_beam[p][s+offset][ch][pol] = UCMPLX4_TO_CMPLX_FLT(data[D_IDX(s,ch,i,nchan,ninput)]);
                detected_beam[p][s+offset][ch][pol] = cuCmul(detected_beam[p][s+offset][ch][pol], make_cuDoubleComplex(128.0, 0.0));
            }
        }
        else // beamform (the default mode)
        {
            cu_form_beam( data, sample_rate, gdelays.d_phi, invJi, timestep_idx,
                    npointing, nant, nchan, npol, outpol_coh, invw, &gf,
                    detected_beam, data_buffer_coh, data_buffer_incoh,
                    streams, opts.out_incoh, nchunk );
        }

        // Invert the PFB, if requested
        if (opts.out_vdif)
        {
            fprintf( log, "[%f] [%lu/%lu] Inverting the PFB (full)\n",
                            NOW-begintime, timestep_idx+1, ntimesteps);
            cu_invert_pfb_ord( detected_beam, timestep_idx, npointing,
                               sample_rate, nchan, npol, vf->sizeof_buffer,
                               &gi, data_buffer_vdif );
        }
        calc_time[timestep_idx] = clock() - start;


        // Write out for each pointing
        for ( p = 0; p < npointing; p++)
        {
            start = clock();
            fprintf( log, "[%f] [%lu/%lu] [%d/%d] Writing data to file(s)\n",
                    NOW-begintime, timestep_idx+1, ntimesteps, p+1, npointing );

            if (opts.out_coh)
                psrfits_write_second( &pf[p], data_buffer_coh, nchan,
                                      outpol_coh, p );
            if (opts.out_incoh && p == 0)
                psrfits_write_second( &pf_incoh[p], data_buffer_incoh,
                                      nchan, outpol_incoh, p );
            if (opts.out_vdif)
                vdif_write_second( &vf[p], &vhdr,
                                   data_buffer_vdif + p * vf->sizeof_buffer );
            write_time[timestep_idx][p] = clock() - start;
        }
    }

    // Calculate total processing times
    float read_sum = 0, delay_sum = 0, calc_sum = 0, write_sum = 0;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        read_sum  += (float) read_time[timestep_idx];
        delay_sum += (float) delay_time[timestep_idx];
        calc_sum  += (float) calc_time[timestep_idx];
        for ( p = 0; p < npointing; p++)
            write_sum += (float) write_time[timestep_idx][p];
    }
    float read_mean, delay_mean, calc_mean, write_mean;
    read_mean  = read_sum  / ntimesteps;
    delay_mean = delay_sum / ntimesteps / npointing;
    calc_mean  = calc_sum  / ntimesteps / npointing;
    write_mean = write_sum / ntimesteps / npointing;

    // Calculate the standard deviations
    float read_std = 0, delay_std = 0, calc_std = 0, write_std = 0;
    for (timestep_idx = 0; timestep_idx < ntimesteps; timestep_idx++)
    {
        read_std  += pow((float)read_time[timestep_idx]  - read_mean,  2);
        delay_std += pow((float)delay_time[timestep_idx] - delay_mean / npointing, 2);
        calc_std  += pow((float)calc_time[timestep_idx]  - calc_mean / npointing,  2);
        for ( p = 0; p < npointing; p++)
            write_std += pow((float)write_time[timestep_idx][p] - write_mean / npointing, 2);
    }
    read_std  = sqrt( read_std  / ntimesteps );
    delay_std = sqrt( delay_std / ntimesteps / npointing );
    calc_std  = sqrt( calc_std  / ntimesteps / npointing );
    write_std = sqrt( write_std / ntimesteps / npointing );

    fprintf( log, "\n[%f]  *****FINISHED BEAMFORMING*****\n\n", NOW-begintime );

    fprintf( log, "[%f]  Total read  processing time: %9.3f s\n",
                     NOW-begintime, read_sum / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Mean  read  processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, read_mean / CLOCKS_PER_SEC, read_std / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Total delay processing time: %9.3f s\n",
                     NOW-begintime, delay_sum / CLOCKS_PER_SEC);
    fprintf( log, "[%f]  Mean  delay processing time: %9.3f +\\- %8.3f s\n",
                     NOW-begintime, delay_mean / CLOCKS_PER_SEC, delay_std / CLOCKS_PER_SEC);
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
    destroy_invJi( invJi, nant, nchan, npol );
    destroy_detected_beam( detected_beam, npointing, 2*sample_rate, nchan );

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

    free( opts.pointings_file );
    free( opts.datadir        );
    free( cal.caldir          );
    free( opts.metafits       );
    free( opts.synth_filter   );

    if (opts.out_incoh)
    {
        free( pf_incoh[0].sub.data        );
        free( pf_incoh[0].sub.dat_freqs   );
        free( pf_incoh[0].sub.dat_weights );
        free( pf_incoh[0].sub.dat_offsets );
        free( pf_incoh[0].sub.dat_scales  );
    }
    for (p = 0; p < npointing; p++)
    {
        if (opts.out_coh)
        {
            free( pf[p].sub.data        );
            free( pf[p].sub.dat_freqs   );
            free( pf[p].sub.dat_weights );
            free( pf[p].sub.dat_offsets );
            free( pf[p].sub.dat_scales  );
        }
        if (opts.out_vdif)
        {
            free( vf[p].b_scales  );
            free( vf[p].b_offsets );
        }
    }
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
    if (cal_metadata != NULL)
        mwalib_metafits_metadata_free( cal_metadata );

    // Clean up memory used for calibration solutions
    free_rts( D, cal_metadata );

    return EXIT_SUCCESS;
}


void usage() {
    printf( "\nusage: make_mwa_tied_array_beam [OPTIONS]\n");

    printf( "\nREQUIRED OPTIONS\n\n"
            "\t-b, --begin=GPSTIME       Begin time of observation, in GPS seconds\n"
            "\t-e, --end=GPSTIME         End time of observation, in GPS seconds\n"
            "\t-P, --pointings=FILE      FILE containing RA and Decs of multiple pointings\n"
            "\t                          in the format hh:mm:ss.s dd:mm:ss.s ...\n"
            "\t-d, --data-location=PATH  PATH is the directory containing the recombined data\n"
            "\t-m, --metafits=FILE       FILE is the metafits file for the target observation\n"
            "\t-f, --coarse-chan=N       Absolute coarse channel number (0-255)\n" );

    printf( "\nOUTPUT OPTIONS\n\n"
            "\t-i, --incoh                Turn on incoherent PSRFITS beam output.                             [default: OFF]\n"
            "\t-p, --psrfits              Turn on coherent PSRFITS output (will be turned on if none of\n"
            "\t                           -i, -p, -u, -v are chosen).                                         [default: OFF]\n"
            "\t-v, --vdif                 Turn on VDIF output with upsampling                                 [default: OFF]\n"
            "\t-s, --summed               Turn on summed polarisations of the coherent output (only Stokes I) [default: OFF]\n"
            "\t-t, --max_t                Maximum number of seconds per output fits file. [default: 200]\n"
            "\t-A, --antpol=ant           Do not beamform. Instead, only operate on the specified ant\n"
            "\t                           stream (0-127)\n" 
            "\t-S, --synth_filter=filter  Apply the named filter during high-time resolution synthesis.    [default: LSQ12]\n"
            "\t                           filter can be MIRROR or LSQ12.\n" );

    printf( "\nCALIBRATION OPTIONS (RTS)\n\n"
            "\t-c, --cal-metafits=FILE  FILE is the metafits file pertaining to the calibration solution\n"
            "\t-C, --cal-location=PATH  PATH is the directory (RTS) or the file (OFFRINGA) containing the calibration solution\n"
            "\t-R, --ref-ant=ANT         Rotate the phases of the XX and YY elements of the calibration\n"
            "\t                          Jones matrices so that the phases of tile ANT align. If ANT is\n"
            "\t                          outside the range 0-127, no phase rotation is done               [default: 0]\n"
            "\t-X, --cross-terms         Retain the XY and YX terms of the calibration solution           [default: off]\n"
            "\t-U, --no-XY-phase         Do not apply the XY phase correction to the calibration solution\n" );

    printf( "\nCALIBRATION OPTIONS (OFFRINGA) -- NOT YET SUPPORTED\n\n"
            "\t-O, --offringa            The calibration solution is in the Offringa format instead of\n"
            "\t                          the default RTS format. In this case, the argument to -C should\n" 
            "\t                          be the full path to the binary solution file.\n" );

    printf( "\nOTHER OPTIONS\n\n"
            "\t-h, --help                Print this help and exit\n"
            "\t-g, --gpu-mem=N     The maximum amount of GPU memory you want make_beam to use in GB [default: -1]\n"
            "\t-V, --version             Print version number and exit\n\n");
}



void make_beam_parse_cmdline(
        int argc, char **argv, struct make_beam_opts *opts, struct calibration *cal )
{
    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"begin",           required_argument, 0, 'b'},
                {"end",             required_argument, 0, 'e'},
                {"incoh",           no_argument,       0, 'i'},
                {"psrfits",         no_argument,       0, 'p'},
                {"vdif",            no_argument,       0, 'v'},
                {"summed",          no_argument,       0, 's'},
                {"max_t",           required_argument, 0, 't'},
                {"synth_filter",    required_argument, 0, 'S'},
                {"antpol",          required_argument, 0, 'A'},
                {"pointings",       required_argument, 0, 'P'},
                {"data-location",   required_argument, 0, 'd'},
                {"cal-location",    required_argument, 0, 'C'},
                {"metafits",        required_argument, 0, 'm'},
                {"cal-metafits",    required_argument, 0, 'c'},
                {"coarse-chan",     required_argument, 0, 'f'},
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
                             "A:b:c:d:C:e:f:g:him:OpP:R:sS:t:UvVX",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'A':
                    opts->out_bf = 0; // Turn off normal beamforming
                    opts->out_ant = atoi(optarg); // 0-127
                    break;
                case 'b':
                    opts->begin = atol(optarg);
                    break;
                case 'c':
                    cal->metafits = strdup(optarg);
                    break;
                case 'C':
                    cal->caldir = strdup(optarg);
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
                    opts->pointings_file = strdup(optarg);
                    break;
                case 'R':
                    cal->ref_ant = atoi(optarg);
                    break;
                case 'S':
                    opts->synth_filter = strdup(optarg);
                    break;
                case 's':
                    opts->out_summed = 1;
                    break;
                case 't':
                    opts->max_sec_per_file = atoi(optarg);
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
                    fprintf( stderr, "error: make_beam_parse_cmdline: "
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
    assert( opts->pointings_file != NULL );
    assert( opts->datadir      != NULL );
    assert( cal->caldir       != NULL );
    assert( opts->metafits     != NULL );
    if (opts->out_coh || opts->out_vdif)
        assert( cal->metafits );

    // If a calibration metafits was supplied, then the default
    // calibration format is RTS
    if (cal->cal_type == CAL_NONE && cal->metafits != NULL)
        cal->cal_type = CAL_RTS;

    // If neither -i, -p, nor -v were chosen, set -p by default
    if ( !opts->out_incoh && !opts->out_coh && !opts->out_vdif )
    {
        opts->out_coh = 1;
    }

    // If the reference antenna is outside the range of antennas, issue
    // a warning and turn off phase rotation.
    /*if (cal->ref_ant < 0 || cal->ref_ant >= opts->nant)
    {
        fprintf( stderr, "warning: tile %d outside of range 0-%d. "
                "Calibration phase rotation turned off.\n",
                cal->ref_ant, opts->nant-1 );
        cal->ref_ant = -1; // This defines the meaning of -1
                                // = don't do phase rotation
    }*/
}



cuDoubleComplex ****create_invJi( int nant, int nchan, int npol )
// Allocate memory for (inverse) Jones matrices
{
    int ant, pol, ch; // Loop variables
    cuDoubleComplex ****invJi;
    invJi = (cuDoubleComplex ****)malloc( nant * sizeof(cuDoubleComplex ***) );

    for (ant = 0; ant < nant; ant++)
    {
        invJi[ant] =(cuDoubleComplex ***)malloc( nchan * sizeof(cuDoubleComplex **) );

        for (ch = 0; ch < nchan; ch++)
        {
            invJi[ant][ch] = (cuDoubleComplex **)malloc( npol * sizeof(cuDoubleComplex *) );

            for (pol = 0; pol < npol; pol++)
                invJi[ant][ch][pol] = (cuDoubleComplex *)malloc( npol * sizeof(cuDoubleComplex) );
        }
    }
    return invJi;
}


void destroy_invJi( cuDoubleComplex ****array, int nant, int nchan, int npol )
{
    int ant, ch, pol;
    for (ant = 0; ant < nant; ant++)
    {
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
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
