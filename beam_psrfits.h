/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef BEAM_PSRFITS_H
#define BEAM_PSRFITS_H

#include <mpi.h>
#include <mwalib.h>
#include "geometry.h"
#include "psrfits.h"


typedef struct mpi_psrfits_t
{
    MPI_Datatype    coarse_chan_spectrum;
    MPI_Datatype    total_spectrum_type;
    MPI_Datatype    spliced_type;

    MPI_Request     request_data;
    MPI_Request     request_offsets;
    MPI_Request     request_scales;

    int             ncoarse_chans;

    struct psrfits  coarse_chan_pf;
    struct psrfits  spliced_pf;

    int             writer_id;
} mpi_psrfits;

void populate_psrfits_header(
        struct psrfits   *pf,
        MetafitsMetadata *obs_metadata,
        VoltageMetadata  *vcs_metadata,
        int               coarse_chan_idx,
        int               max_sec_per_file,
        int               outpol,
        struct beam_geom *beam_geom_vals,
        char             *incoh_basename,
        bool              is_coherent );

void populate_spliced_psrfits_header(
        struct psrfits   *pf,
        MetafitsMetadata *obs_metadata,
        VoltageMetadata  *vcs_metadata,
        int               first_coarse_chan_idx,
        int               ncoarse_chans,
        int               max_sec_per_file,
        int               outpol,
        struct beam_geom *beam_geom_vals,
        char             *incoh_basename,
        bool              is_coherent );

void free_psrfits( struct psrfits *pf );

float *create_data_buffer_psrfits( size_t size );

void init_mpi_psrfits(
        mpi_psrfits *mpf,
        MetafitsMetadata *obs_metadata,
        VoltageMetadata *vcs_metadata,
        int world_size,
        int world_rank,
        int max_sec_per_file,
        int nstokes,
        struct beam_geom *bg,
        char *outfile,
        int writer_id,
        bool is_coherent );

void free_mpi_psrfits( mpi_psrfits *mpf );

void gather_splice_psrfits( mpi_psrfits *mpf );
void wait_splice_psrfits( mpi_psrfits *mpf );

#endif
