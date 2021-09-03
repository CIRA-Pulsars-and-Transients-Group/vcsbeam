/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __METADATA_H__
#define __METADATA_H__

#define ERROR_MESSAGE_LEN  1024
#define MAX_COMMAND_LENGTH 1024

#include <mwalib.h>

/* The following is a patch to provide constants which are not available in mwalib */
#ifndef SPEED_OF_LIGHT_IN_VACUUM_M_PER_S
#define SPEED_OF_LIGHT_IN_VACUUM_M_PER_S 299792458.0
#endif

#ifndef MWA_LATITUDE_RADIANS
#define MWA_LATITUDE_RADIANS -0.4660608448386394
#endif

#ifndef MWA_LONGITUDE_RADIANS
#define MWA_LONGITUDE_RADIANS 2.0362898668561042
#endif

#ifndef MWA_ALTITUDE_METRES
#define MWA_ALTITUDE_METRES 377.827
#endif
/* End replacement constants */

char **create_filenames(
        const struct MetafitsContext  *metafits_context,
        const struct MetafitsMetadata *metafits_metadata,
        unsigned long int              begin_gps,
        unsigned long int              nseconds,
        char                          *datadir,
        uintptr_t                      begin_coarse_chan_idx,
        uintptr_t                      ncoarse_chans
        );

void destroy_filenames( char **filenames, int nfiles );

void get_mwalib_metafits_metadata(
        char              *filename,
        MetafitsMetadata **metadata,
        MetafitsContext  **context
        );

void get_mwalib_voltage_metadata(
        VoltageMetadata  **vcs_metadata,
        VoltageContext   **vcs_context,
        MetafitsMetadata **obs_metadata,
        MetafitsContext   *obs_context,
        unsigned long int  begin_gps,
        int                nseconds,
        char               *datadir,
        uintptr_t          coarse_chan_idx,
        int               *ncoarse_chans
        );

long unsigned int get_relative_gps( MetafitsMetadata *obs_metadata, long int relative_begin );
long unsigned int parse_begin_string( MetafitsMetadata *obs_metadata, char *begin_str );
uintptr_t parse_coarse_chan_string( MetafitsMetadata *obs_metadata, char *begin_coarse_chan_str );

#endif
