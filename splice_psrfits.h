/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __SPLICE_PSRFITS_H__
#define __SPLICE_PSRFITS_H__

void splice_psrfits_usage();
int prep_data( void *to_load );
int cleanup( void *to_clean );
float round_three( float var );
int append( char **to_append, void *total, int n );
int splice_psrfits_main( int argc, char *argv[] );

#endif
