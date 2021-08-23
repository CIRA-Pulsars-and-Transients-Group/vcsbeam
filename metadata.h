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

char **create_filenames(
        const struct MetafitsContext  *metafits_context,
        const struct MetafitsMetadata *metafits_metadata,
        unsigned long int              begin,
        unsigned long int              nseconds,
        char                          *datadir,
        uintptr_t                      rec_channel
        );

void destroy_filenames( char **filenames, int nfiles );

void get_mwalib_metafits_metadata(
        char              *filename,
        MetafitsMetadata **metadata,
        MetafitsContext  **context
        );

void get_mwalib_metadata(
        MetafitsMetadata **obs_metadata,
        VoltageMetadata  **vcs_metadata,
        VoltageContext   **vcs_context,
        char              *obs_metafits,
        unsigned long int  begin_gps,
        int                nseconds,
        char               *datadir,
        uintptr_t          rec_channel
        );

long unsigned int get_relative_gps( MetafitsMetadata *obs_metadata, long int relative_begin );
long unsigned int parse_begin_string( MetafitsMetadata *obs_metadata, char *begin_str );

#endif
