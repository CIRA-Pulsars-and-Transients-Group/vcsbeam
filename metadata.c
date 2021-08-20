/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mwalib.h>
#include "make_beam.h"
#include "metadata.h"

char **create_filenames(
        const struct MetafitsContext  *metafits_context,
        const struct MetafitsMetadata *metafits_metadata,
        unsigned long int              begin,
        unsigned long int              end,
        char                          *datadir,
        uintptr_t                      rec_channel
        )
// Create an array of filenames; free with destroy_filenames()
{
    // Buffer for mwalib error messages
    char error_message[ERROR_MESSAGE_LEN];

    // Calculate the number of files
    unsigned int ntimesteps = end - begin + 1;
    if (ntimesteps <= 0) {
        fprintf( stderr, "Cannot beamform on %d files (between %ld and %ld)\n",
                 ntimesteps, begin, end);
        exit(EXIT_FAILURE);
    }
    // Allocate memory for the file name list
    char filename[MAX_COMMAND_LENGTH]; // Just the mwalib-generated filename (without the path)
    char **filenames = (char **)malloc( ntimesteps*sizeof(char *) ); // The full array of filenames, including the paths

    // Get the coarse channel index
    unsigned int coarse_chan_idx;
    for (coarse_chan_idx = 0; coarse_chan_idx < metafits_metadata->num_metafits_coarse_chans; coarse_chan_idx++)
        if (metafits_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number == (uintptr_t)rec_channel)
            break;

    // Allocate memory and write filenames
    unsigned int timestep, timestep_idx, filename_idx;
    for (timestep = begin; timestep <= end; timestep++)
    {
        timestep_idx = timestep - (metafits_metadata->metafits_timesteps[0].gps_time_ms / 1000); // <-- This assumes timesteps are always contiguous
        filename_idx = timestep - begin;

        filenames[filename_idx] = (char *)malloc( MAX_COMMAND_LENGTH*sizeof(char) );
        if (mwalib_metafits_get_expected_volt_filename(
                    metafits_context,
                    timestep_idx,
                    coarse_chan_idx,
                    filename,
                    MAX_COMMAND_LENGTH,
                    error_message,
                    ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
        {
            fprintf( stderr, "error (mwalib): cannot create filenames: %s\n", error_message );
            exit(EXIT_FAILURE);
        }
        sprintf( filenames[filename_idx], "%s/%s",
                 datadir, filename );
    }

    return filenames;
}

void destroy_filenames( char **filenames, int nfiles )
{
    int i;
    for (i = 0; i < nfiles; i++)
        free( filenames[i] );
    free( filenames );
}


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
        )
/* Create the metadata structs using mwalib's API.
 *
 * obs_metadata applies to the metafits metadata for the target observation, and cannot be null
 * vcs_metadata applies to the voltage metadata for the target observation, and cannot be null
 * cal_metadata applies to the metafits metadata for the calibration observation, and may be NULL
 *     (if only an incoherent beam is requested)
 * vcs_context applies to the voltage context for the target observation, and cannot be null
 */
{
    char error_message[ERROR_MESSAGE_LEN];
    error_message[0] = '\0'; // <-- Just to avoid a compiler warning about uninitialised variables

    // Each metadata is constructed from mwalib "contexts"
    // These are only temporarily needed to construct the other metadata/context structs,
    // and will be freed at the end of this function
    MetafitsContext *obs_context = NULL;
    MetafitsContext *cal_context = NULL;

    // First, get the metafits context for the given observation metafits file in order to create
    // a list of filenames for creating the voltage context. This metafits context will be remade
    // from the voltage context later, in order to get the antenna ordering correct (which depends
    // on the voltage type)

    // Create OBS_CONTEXT
    if (mwalib_metafits_context_new2( obs_metafits, &obs_context, error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( obs_context, NULL, NULL, obs_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create list of filenames
    char **filenames = create_filenames( obs_context, *obs_metadata, begin, end, datadir, rec_channel );

    // Now that we have the filenames, we don't need this version of the obs context and metadata.
    mwalib_metafits_metadata_free( *obs_metadata );
    mwalib_metafits_context_free( obs_context );

    // Create an mwalib voltage context, voltage metadata, and new obs metadata (now with correct antenna ordering)
    // (MWALIB is expecting a const array, so we will give it one!)
    unsigned int ntimesteps = end - begin + 1;
    const char **voltage_files = (const char **)malloc( sizeof(char *) * ntimesteps );
    unsigned int i;
    for (i = 0; i < ntimesteps; i++)
        voltage_files[i] = filenames[i];

    // Create VCS_CONTEXT
    if (mwalib_voltage_context_new( obs_metafits, voltage_files, ntimesteps, vcs_context, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( NULL, NULL, *vcs_context, obs_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata from voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create VCS_METADATA
    if (mwalib_voltage_metadata_get( *vcs_context, vcs_metadata, error_message, ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot get metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Free memory
    destroy_filenames( filenames, ntimesteps );
    free( voltage_files );

    // Can stop at this point if no calibration metadata is requested
    if (cal_metadata == NULL)
        return;

    // Create CAL_CONTEXT
    if (mwalib_metafits_context_new2( cal_metafits, &cal_context, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create cal metafits context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create CAL_METADATA
    if (mwalib_metafits_metadata_get( cal_context, NULL, NULL, cal_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create cal metafits metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }
}
