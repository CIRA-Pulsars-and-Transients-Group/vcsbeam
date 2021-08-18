/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __PRIMARY_BEAM_H__
#define __PRIMARY_BEAM_H__

#include <mwa_hyperbeam.h>
#include <mwalib.h>


void create_delays_amps_from_metafits( MetafitsMetadata *metafits_metadata, uint32_t ***delays, double ***amps );
void free_delays_amps( MetafitsMetadata *metafits_metadata, uint32_t **delays, double **amps );

void parallactic_angle_correction(
    double *P,    // output rotation matrix
    double lat,   // observing latitude (radians)
    double az,    // azimuth angle (radians)
    double za);   // zenith angle (radians)

int hash_dipole_config( double * );

#endif
