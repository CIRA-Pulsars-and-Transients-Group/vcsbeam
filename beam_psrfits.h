/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef BEAM_PSRFITS_H
#define BEAM_PSRFITS_H

#include <mwalib.h>
#include "psrfits.h"

void printf_psrfits( struct psrfits *pf );  /* Prints values in psrfits struct to stdout */

void populate_psrfits_header(
        MetafitsMetadata *obs_metadata,
        VoltageMetadata  *vcs_metadata,
        int               coarse_chan_idx,
        struct psrfits    pf[],
        int               max_sec_per_file,
        int               outpol,
        struct beam_geom *beam_geom_vals,
        int               npointing,
        bool              is_coherent );

void correct_psrfits_stt( struct psrfits *pf );

void psrfits_write_second( struct psrfits *pf, float *data_buffer, int nchan,
        int outpol, int p);


#endif
