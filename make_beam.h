/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef MAKE_BEAM_SMALL_H
#define MAKE_BEAM_SMALL_H

#include <stdlib.h>
#include "beam_common.h"
#include <cuComplex.h>


#define MAX_COMMAND_LENGTH 1024

void usage();
void make_beam_parse_cmdline( int argc, char **argv, struct make_beam_opts *opts );

char **create_filenames( struct make_beam_opts *opts, int obsid );
void  destroy_filenames( char **filenames, struct make_beam_opts *opts );

cuDoubleComplex ****create_complex_weights( int npointing, int nstation, int nchan, int npol );
void              destroy_complex_weights( cuDoubleComplex ****array, int npointing,
                                           int nstation, int nchan );

cuDoubleComplex ****create_invJi( int nstation, int nchan, int pol );
void              destroy_invJi( cuDoubleComplex ****array, int nstation, int nchan, int npol );

cuDoubleComplex ****create_detected_beam( int npointing, int nsamples, int nchan, int npol );
void              destroy_detected_beam( cuDoubleComplex ****array, int npointing,
                                         int nsamples, int nchan );

float *create_data_buffer_psrfits( size_t size );
float *create_data_buffer_vdif( size_t size );

#endif
