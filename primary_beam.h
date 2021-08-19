/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __PRIMARY_BEAM_H__
#define __PRIMARY_BEAM_H__

#include <mwa_hyperbeam.h>
#include <mwalib.h>
#include <cuComplex.h>
#include "beam_common.h"

#define NCONFIGS 138

void create_delays_amps_from_metafits(
        MetafitsMetadata *metafits_metadata, uint32_t ***delays, double ***amps );

void free_delays_amps(
        MetafitsMetadata *metafits_metadata, uint32_t **delays, double **amps );

void parallactic_angle_correction(
    double *P,    // output rotation matrix
    double lat,   // observing latitude (radians)
    double az,    // azimuth angle (radians)
    double za);   // zenith angle (radians)

int hash_dipole_config( double * );

void calc_primary_beam(
        cuDoubleComplex ***B,
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        int                coarse_chan_idx,
        struct beam_geom  *beam_geom_vals,
        FEEBeam           *beam,
        uint32_t         **delays,
        double           **amps,
        uintptr_t          npointings );

cuDoubleComplex ***malloc_primary_beam(
        MetafitsMetadata  *obs_metadata,
        uintptr_t          npointings );

void free_primary_beam(
        cuDoubleComplex ***B,
        MetafitsMetadata  *obs_metadata,
        uintptr_t          npointings );

#endif
