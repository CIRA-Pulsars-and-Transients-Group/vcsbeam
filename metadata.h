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
        unsigned long int              end,
        char                          *datadir,
        uintptr_t                      rec_channel
        );

void destroy_filenames( char **filenames, int nfiles );

void get_mwalib_metadata(
        MetafitsMetadata **obs_metadata,
        VoltageMetadata  **vcs_metadata,
        VoltageContext   **vcs_context,
        MetafitsMetadata **cal_metadata,
        char              *obs_metafits,
        char              *cal_metafits,
        unsigned long int  begin,
        unsigned long int  end,
        char               *datadir,
        uintptr_t          rec_channel
        );

#endif
