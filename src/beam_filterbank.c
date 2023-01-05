/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <sigproc/sigproc.h>

int sigproc_test_function()
{
    printf( "here! little endian? %d\n", little_endian() );
}

/*
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
*/

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
/*
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
*/

/**
 * Write VDIF data to file.
 *
 * @param vf A struct containing metadata about the target VDIF file
 * @param output The buffer to be written to file
 *
 * @todo Make sure that the header file is unneccesarily re-written each time
 * this function is called. If so, put this bit somewhere more sensible.
 */
/*
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
*/


