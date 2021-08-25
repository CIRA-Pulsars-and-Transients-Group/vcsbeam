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
#include "jones.h"

#define NCONFIGS           138
#define DEAD_CONFIG        (NCONFIGS - 1)
#define MANY_DEAD_DIPOLES  -1

#define PB_IDX(p,a,pol,na,npol) ((p)*(na)*(npol) + (a)*(npol) + (pol))

typedef struct primary_beam_t
{
    cuDoubleComplex  *B;
    FEEBeam          *beam;
    uint32_t        **delays;
    double          **amps;
    uintptr_t         npointings;
    uintptr_t         nant;
    uintptr_t         npol;
    uint32_t          freq_hz;
    MetafitsMetadata *obs_metadata;
} primary_beam;

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
        primary_beam      *pb,
        struct beam_geom  *beam_geom_vals );

void create_primary_beam(
        primary_beam      *pb,
        MetafitsMetadata  *obs_metadata,
        uintptr_t          coarse_chan,
        uintptr_t          npointings );

void free_primary_beam( primary_beam *pb );

#endif
