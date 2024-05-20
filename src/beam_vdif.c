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
#include <vdifio.h>

#include "vcsbeam.h"

#include "mwa_header.h"
#include "ascii_header.h"

/**
 * Converts floats to 8-bit integers.
 *
 * @param[in]  f   The (source) buffer of floats
 * @param      n   The number of floats in `f`
 * @param      min The smallest allowed float
 * @param      max The largest allowed float
 * @param[out] i   The (destination) buffer of 8-bit integers
 */
void float2int8_trunc(float *f, int n, float min, float max, int8_t *i)
{
    int j;
    for (j = 0; j < n; j++) {
        f[j] = (f[j] > max) ? (max) : f[j];
        f[j] = (f[j] < min) ? (min) : f[j];
        i[j] = (int8_t) rint(f[j]);

    }
}

/**
 * Writes a second's worth of data to a VDIF file
 *
 * @param vf A struct containing metadata about the target VDIF file
 * @param vhdr A pointer to the VDIF frame header buffer
 * @param data_buffer_vdif The data to be written out
 *
 * The data to be written out are first normalised and scaled to fit
 * in the range -126 to 127 and demoted to integers.
 * "Blocks" of data are then written to file, along with the appropriate
 * binary headers for each "frame", via the function vdif_write_data().
 */
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
    /*
    for (int i = 0; i < 96; i++)
    {
        printf("%02x ", *(uint8_t *)(out_buffer_8_vdif + i));
        if ((i % 16) == 15) printf("\n");
        if ((i % 16) == 7) printf(" ");
    }
    printf("...\n");
    */

    free( out_buffer_8_vdif );
}

/**
 * Write VDIF data to file.
 *
 * @param vf A struct containing metadata about the target VDIF file
 * @param output The buffer to be written to file
 *
 * @todo Make sure that the header file is unneccesarily re-written each time
 * this function is called. If so, put this bit somewhere more sensible.
 */
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
    
    ascii_header_set( ascii_header, "TELESCOPE",  "%s", vf->telescope   );
    ascii_header_set( ascii_header, "MODE",       "%s", vf->obs_mode    );
    ascii_header_set( ascii_header, "INSTRUMENT", "%s", "VDIF"          );
    ascii_header_set( ascii_header, "DATAFILE",   "%s", filename        );

    ascii_header_set( ascii_header, "UTC_START",  "%s", vf->date_obs    );
    ascii_header_set( ascii_header, "MJD_START",  "%f", vf->MJD_start   );

    ascii_header_set( ascii_header, "SEC_OFFSET", "%f", vf->sec_offset  );
    ascii_header_set( ascii_header, "MJD_EPOCH",  "%f", vf->MJD_epoch   );

    ascii_header_set( ascii_header, "SOURCE",     "%s", vf->source      );
    ascii_header_set( ascii_header, "RA",         "%s", vf->ra_str      );
    ascii_header_set( ascii_header, "DEC",        "%s", vf->dec_str     );

    ascii_header_set( ascii_header, "FREQ",       "%f", vf->fctr        );
    ascii_header_set( ascii_header, "BW",         "%f", vf->BW          );
    ascii_header_set( ascii_header, "TSAMP",      "%f", 1e6/(float) vf->sample_rate );

    ascii_header_set( ascii_header, "NBIT",       "%d", vf->bits        );
    ascii_header_set( ascii_header, "NDIM",       "%d", vf->iscomplex+1 );
    ascii_header_set( ascii_header, "NPOL",       "%d", vf->npol        );

    sprintf( filename, "%s.hdr", vf->basefilename );
    fs = fopen( filename,"w" );
    fwrite( ascii_header, MWA_HEADER_SIZE, 1, fs );
    fclose( fs );

}

/**
 * Populates a VDIF header with data derived from the observation.
 *
 * @param vm The VCSBeam context struct
 * @param beam_geom_vals A `beam_geom` struct containing pointing information
 * @param mjd_start The start MJD of the observation
 * @param sec_offset The offset from the start of the observation in seconds
 */
void vmPopulateVDIFHeader(
        vcsbeam_context  *vm,
        beam_geom        *beam_geom_vals,
        double           mjd_start,
        double           sec_offset )
{
    // Write log message
    sprintf( vm->log_message, "Preparing headers for output (receiver channel %lu)",
            vm->obs_metadata->metafits_coarse_chans[vm->coarse_chan_idxs_to_process[0]].rec_chan_number );
    logger_timed_message( vm->log, vm->log_message );

    vm->vf = (struct vdifinfo *)malloc(vm->npointing * sizeof(struct vdifinfo));

    // Shorthand variables
    int coarse_chan_idx = vm->coarse_chan_idxs_to_process[0];

    // Convert the UTC obs time into a string
    struct tm *ts = gmtime( &(vm->obs_metadata->sched_start_utc) );
    char   time_utc[24];
    strftime( time_utc, sizeof(time_utc), "%Y-%m-%dT%H:%M:%S", ts );

    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        // Define DataFrame dimensions
        vm->vf[p].bits              = 8;   // this is because it is all the downstream apps support (dspsr/diFX)
        vm->vf[p].iscomplex         = 1;   // (it is complex data)
        vm->vf[p].nchan             = 2;   // I am hardcoding this to 2 channels per thread - one per pol
        vm->vf[p].samples_per_frame = 128; // Hardcoding to 128 time-samples per frame
        vm->vf[p].sample_rate       = vm->fine_sample_rate * vm->nfine_chan; // For coarse channelised data
        vm->vf[p].BW                = vm->obs_metadata->coarse_chan_width_hz / 1e6;  // (MHz)

        vm->vf[p].dataarraylength = (vm->vf[p].nchan * (vm->vf[p].iscomplex+1) * vm->vf[p].samples_per_frame); // = 512
        vm->vf[p].frame_length  = vm->vf[p].dataarraylength + VDIF_HEADER_BYTES;                       // = 544
        vm->vf[p].threadid      = 0;
        sprintf( vm->vf[p].stationid, "mw" );

        vm->vf[p].frame_rate = vm->fine_sample_rate;                                                      // = 10000
        vm->vf[p].block_size = vm->vf[p].frame_length * vm->vf[p].frame_rate;                              // = 5440000

        // A single frame (128 samples). Remember vf.nchan is kludged to npol
        vm->vf[p].sizeof_beam = vm->vf[p].samples_per_frame * vm->vf[p].nchan * (vm->vf[p].iscomplex+1);       // = 512

        // One full second (1.28 million 2 bit samples)
        vm->vf[p].sizeof_buffer = vm->vf[p].frame_rate * vm->vf[p].sizeof_beam;                            // = 5120000

        createVDIFHeader( &vm->vhdr, vm->vf[p].dataarraylength, vm->vf[p].threadid, vm->vf[p].bits, vm->vf[p].nchan,
                                vm->vf[p].iscomplex, vm->vf[p].stationid);

        // Now we have to add the time
        uint64_t start_day = beam_geom_vals->intmjd;
        uint64_t start_sec = roundf( beam_geom_vals->fracmjd * 86400.0 );
        uint64_t mjdsec    = (start_day * 86400) + start_sec; // Note the VDIFEpoch is strange - from the standard

        setVDIFEpochMJD( &vm->vhdr, start_day );
        setVDIFFrameMJDSec( &vm->vhdr, mjdsec );
        setVDIFFrameNumber( &vm->vhdr, 0 );

        strcpy( vm->vf[p].exp_name, vm->obs_metadata->project_id );
        snprintf( vm->vf[p].scan_name, 17, "%d", vm->obs_metadata->obs_id );

        vm->vf[p].b_scales   = (float *)malloc( sizeof(float) * vm->vf[p].nchan );
        vm->vf[p].b_offsets  = (float *)malloc( sizeof(float) * vm->vf[p].nchan );
        vm->vf[p].got_scales = 1;

        strncpy( vm->vf[p].telescope, "MWA", 24);
        strncpy( vm->vf[p].obs_mode,  "PSR", 8);

        // Determine the RA and Dec strings
        double ra2000  = beam_geom_vals[p].mean_ra  * R2D;
        double dec2000 = beam_geom_vals[p].mean_dec * R2D;

        dec2hms(vm->vf[p].ra_str,  ra2000/15.0, 0); // 0 = no '+' sign
        dec2hms(vm->vf[p].dec_str, dec2000,     1); // 1 = with '+' sign

        strncpy( vm->vf[p].date_obs, time_utc, sizeof(time_utc) );

        vm->vf[p].MJD_start  = mjd_start;
        vm->vf[p].sec_offset = sec_offset;
        vm->vf[p].MJD_epoch  = beam_geom_vals->intmjd + beam_geom_vals->fracmjd;
        vm->vf[p].fctr       = vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_centre_hz / 1e6; // (MHz)
        vm->vf[p].npol       = 2;
        strncpy( vm->vf[p].source, "unset", 24 );

        // The output file basename
        sprintf( vm->vf[p].basefilename, "%s_%s_%s_%s_ch%03ld",
                 vm->vf[p].exp_name,
                 vm->vf[p].scan_name,
                 vm->vf[p].ra_str,
                 vm->vf[p].dec_str,
                 vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
    }
}

/**
 * Convert from two's complement to "offset binary" (??)
 *
 * @param i An array of integers to be converted
 * @param n The number of integers in `i`
 */
void to_offset_binary(int8_t *i, int n)
{
    int j;
    for (j = 0; j < n; j++) {
        i[j] = i[j] ^ 0x80;
    }
}

