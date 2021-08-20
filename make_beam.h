/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef MAKE_BEAM_SMALL_H
#define MAKE_BEAM_SMALL_H

#include <stdlib.h>
#include "beam_common.h"
#include "calibration.h"
#include <cuComplex.h>


#define MAX_COMMAND_LENGTH 1024

void usage();
void make_beam_parse_cmdline( int argc, char **argv, struct make_beam_opts *opts, struct calibration *cal );


#endif
