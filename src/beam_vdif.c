/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuComplex.h>

#include <mwalib.h>
#include <psrfits.h>
#include <star/pal.h>
#include <star/palmac.h>
#include <vdifio.h>

#include "vcsbeam.h"

#include "mwa_header.h"
#include "vdifio.h"
#include "ascii_header.h"


void float2int8_trunc(float *f, int n, float min, float max, int8_t *i)
{
    int j;
    for (j = 0; j < n; j++) {
        f[j] = (f[j] > max) ? (max) : f[j];
        f[j] = (f[j] < min) ? (min) : f[j];
        i[j] = (int8_t) rint(f[j]);

    }
}


void vdif_write_second( struct vdifinfo *vf, vdif_header *vhdr,
                        float *data_buffer_vdif )
{
    float *data_buffer_ptr = data_buffer_vdif;
    size_t offset_out_vdif = 0;

    int8_t *out_buffer_8_vdif = (int8_t *)malloc(vf->block_size);

    while  (offset_out_vdif < vf->block_size) {

        // Add the current header
        memcpy( (out_buffer_8_vdif + offset_out_vdif), vhdr, VDIF_HEADER_BYTES );

        // Offset into the output array
        offset_out_vdif += VDIF_HEADER_BYTES;

        // Convert from float to int8
        float2int8_trunc( data_buffer_ptr, vf->sizeof_beam, -126.0, 127.0,
                          (out_buffer_8_vdif + offset_out_vdif) );
        to_offset_binary( (out_buffer_8_vdif + offset_out_vdif),
                          vf->sizeof_beam );

        offset_out_vdif += vf->frame_length - VDIF_HEADER_BYTES; // increment output offset
        data_buffer_ptr += vf->sizeof_beam;
        nextVDIFHeader( vhdr, vf->frame_rate );
    }

    // Write a full second's worth of samples
    vdif_write_data( vf, out_buffer_8_vdif );

    free( out_buffer_8_vdif );
}

void vdif_write_data( struct vdifinfo *vf, int8_t *output )
{
    // form the filename
    // there is a standard naming convention
    char  filename[1030];
    sprintf( filename, "%s.vdif", vf->basefilename );

    //fprintf(stderr,"Attempting to open VDIF file for writing: %s\n",filename);
    FILE *fs = fopen( filename, "a" );
    fwrite( output, vf->block_size, 1, fs );
    fclose( fs );

    // write a CPSR2 test header for DSPSR
    char ascii_header[MWA_HEADER_SIZE] = MWA_HEADER_INIT;
    //ascii_header_set( ascii_header, "UTC_START", "%s", vf->date_obs  );
    ascii_header_set( ascii_header, "DATAFILE",   "%s", filename      );
    ascii_header_set( ascii_header, "INSTRUMENT", "%s", "VDIF"        );
    ascii_header_set( ascii_header, "TELESCOPE",  "%s", vf->telescope );
    ascii_header_set( ascii_header, "MODE",       "%s", vf->obs_mode  );
    ascii_header_set( ascii_header, "FREQ",       "%f", vf->fctr      );

    ascii_header_set( ascii_header, "BW",         "%f", vf->BW        );
    ascii_header_set( ascii_header, "RA",         "%s", vf->ra_str    );
    ascii_header_set( ascii_header, "DEC",        "%s", vf->dec_str   );
    ascii_header_set( ascii_header, "SOURCE",     "%s", vf->source    );

    sprintf( filename, "%s.hdr", vf->basefilename );
    fs = fopen( filename,"w" );
    fwrite( ascii_header, MWA_HEADER_SIZE, 1, fs );
    fclose( fs );

}


void populate_vdif_header(
        struct vdifinfo  *vf,
        vdif_header      *vhdr,
        MetafitsMetadata *obs_metadata,
        VoltageMetadata *vcs_metadata,
        int               coarse_chan_idx,
        struct beam_geom *beam_geom_vals,
        int               npointing )
{
    // Convert the UTC obs time into a string
    struct tm *ts = gmtime( &(obs_metadata->sched_start_utc) );
    char   time_utc[24];
    strftime( time_utc, sizeof(time_utc), "%Y-%m-%dT%H:%M:%S", ts );

    // Get the sample rate
    unsigned int sample_rate = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;

    int p;
    for (p = 0; p < npointing; p++)
    {
        // Define DataFrame dimensions
        vf[p].bits              = 8;   // this is because it is all the downstream apps support (dspsr/diFX)
        vf[p].iscomplex         = 1;   // (it is complex data)
        vf[p].nchan             = 2;   // I am hardcoding this to 2 channels per thread - one per pol
        vf[p].samples_per_frame = 128; // Hardcoding to 128 time-samples per frame
        vf[p].sample_rate       = sample_rate * vcs_metadata->num_fine_chans_per_coarse; // For coarse channelised data
        vf[p].BW                = obs_metadata->coarse_chan_width_hz / 1e6;  // (MHz)

        vf[p].dataarraylength = (vf[p].nchan * (vf[p].iscomplex+1) * vf[p].samples_per_frame); // = 512
        vf[p].frame_length  = vf[p].dataarraylength + VDIF_HEADER_BYTES;                       // = 544
        vf[p].threadid      = 0;
        sprintf( vf[p].stationid, "mw" );

        vf[p].frame_rate = sample_rate;                                                        // = 10000
        vf[p].block_size = vf[p].frame_length * vf[p].frame_rate;                              // = 5440000

        // A single frame (128 samples). Remember vf.nchan is kludged to npol
        vf[p].sizeof_beam = vf[p].samples_per_frame * vf[p].nchan * (vf[p].iscomplex+1);       // = 512

        // One full second (1.28 million 2 bit samples)
        vf[p].sizeof_buffer = vf[p].frame_rate * vf[p].sizeof_beam;                            // = 5120000

        createVDIFHeader( vhdr, vf[p].dataarraylength, vf[p].threadid, vf[p].bits, vf[p].nchan,
                                vf[p].iscomplex, vf[p].stationid);

        // Now we have to add the time
        uint64_t start_day = beam_geom_vals->intmjd;
        uint64_t start_sec = roundf( beam_geom_vals->fracmjd * 86400.0 );
        uint64_t mjdsec    = (start_day * 86400) + start_sec; // Note the VDIFEpoch is strange - from the standard

        setVDIFEpochMJD( vhdr, start_day );
        setVDIFFrameMJDSec( vhdr, mjdsec );
        setVDIFFrameNumber( vhdr, 0 );

        strcpy( vf[p].exp_name, obs_metadata->project_id );
        snprintf( vf[p].scan_name, 17, "%d", obs_metadata->obs_id );

        vf[p].b_scales   = (float *)malloc( sizeof(float) * vf[p].nchan );
        vf[p].b_offsets  = (float *)malloc( sizeof(float) * vf[p].nchan );
        vf[p].got_scales = 1;

        strncpy( vf[p].telescope, "MWA", 24);
        strncpy( vf[p].obs_mode,  "PSR", 8);

        // Determine the RA and Dec strings
        double ra2000  = beam_geom_vals[p].mean_ra  * PAL__DR2D;
        double dec2000 = beam_geom_vals[p].mean_dec * PAL__DR2D;

        dec2hms(vf[p].ra_str,  ra2000/15.0, 0); // 0 = no '+' sign
        dec2hms(vf[p].dec_str, dec2000,     1); // 1 = with '+' sign

        strncpy( vf[p].date_obs, time_utc, sizeof(time_utc) );

        vf[p].MJD_epoch = beam_geom_vals->intmjd + beam_geom_vals->fracmjd;
        vf[p].fctr      = obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_centre_hz / 1e6; // (MHz)
        strncpy( vf[p].source, "unset", 24 );

        // The output file basename
        sprintf( vf[p].basefilename, "%s_%s_%s_%s_ch%03ld",
                 vf[p].exp_name,
                 vf[p].scan_name,
                 vf[p].ra_str,
                 vf[p].dec_str,
                 obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
    }
}


cuFloatComplex get_std_dev_complex(cuFloatComplex *input, int nsamples)
{
    // assume zero mean
    float rtotal = 0;
    float itotal = 0;
    float isigma = 0;
    float rsigma = 0;

    int i;
    for (i = 0; i < nsamples; i++)
    {
         rtotal = rtotal+(cuCrealf(input[i])*cuCrealf(input[i]));
         itotal = itotal+(cuCimagf(input[i])*cuCimagf(input[i]));
    }
    rsigma = sqrtf((1.0/(nsamples-1))*rtotal);
    isigma = sqrtf((1.0/(nsamples-1))*itotal);

    return make_cuFloatComplex( rsigma, isigma );
}

void set_level_occupancy(cuFloatComplex *input, int nsamples, float *new_gain)
{
    //float percentage = 0.0;
    //float occupancy = 17.0;
    //float limit = 0.00001;
    //float step = 0.001;
    int i = 0;
    float gain = *new_gain;

    float percentage_clipped = 100;
    //while (percentage_clipped > 0 && percentage_clipped > limit) {
        int clipped = 0;
        for (i = 0; i < nsamples; i++) {
            if (isnan(cuCrealf(input[i])) || isnan(cuCimagf(input[i])))
            {
                fprintf( stderr, "error: set_level_occupancy: input[%d] = "
                                 "NaN\n", i );
                exit(EXIT_FAILURE);
            }
            if (fabs(gain*cuCrealf(input[i])) > 127 || fabs(gain*cuCimagf(input[i])) > 127 )
            {
                clipped++;
            }
        }
        percentage_clipped = ((float) clipped/nsamples) * 100;
        //The reduction in the gain was commented out until we work our a robust solution
        //if (percentage_clipped > limit) {
        //    gain = gain - step;
        //}
        if (clipped > 0)
        {
            fprintf(stdout,"warning: percentage samples clipped %f percent\n",percentage_clipped);
        }
    //}
    *new_gain = gain;
}


void get_mean_complex( cuFloatComplex *input, int nsamples, float *rmean,
                       float *imean, cuFloatComplex *cmean)
{
    int i;

    float rtotal = 0;
    float itotal = 0 ;

    cuFloatComplex ctotal = make_cuFloatComplex( 0.0, 0.0 );

    for (i = 0; i < nsamples; i++)
    {
//if (isnan(cuCrealf(input[i])) || isnan(cuCimagf(input[i]))) { fprintf(stderr, "\ninput[%d] = %e + %e*I\n\n", i, cuCrealf(input[i]), cuCimagf(input[i])); exit(1); }
        rtotal += cuCrealf( input[i] );
        itotal += cuCimagf( input[i] );
        ctotal  = cuCaddf( ctotal, input[i] );
    }

    *rmean = rtotal / nsamples;
    *imean = itotal / nsamples;
    *cmean = make_cuFloatComplex( cuCrealf(ctotal)/(float)nsamples, cuCimagf(ctotal)/(float)nsamples );
}

void normalise_complex(cuFloatComplex *input, int nsamples, float scale)
{
    int i;
    for (i = 0; i < nsamples; i++)
    {
        input[i] = make_cuFloatComplex( cuCrealf(input[i])*scale, cuCimagf(input[i])*scale );
    }
}


void to_offset_binary(int8_t *i, int n)
{
    int j;
    for (j = 0; j < n; j++) {
        i[j] = i[j] ^ 0x80;
    }
}

float *create_data_buffer_vdif( size_t size )
{
    float *ptr  = (float *)malloc( size * sizeof(float) );
    return ptr;
}

