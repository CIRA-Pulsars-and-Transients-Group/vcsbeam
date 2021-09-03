/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mwalib.h>
#include "metadata.h"

char **create_filenames(
        const struct MetafitsContext  *metafits_context,
        const struct MetafitsMetadata *metafits_metadata,
        unsigned long int              begin_gps,
        unsigned long int              nseconds,
        char                          *datadir,
        uintptr_t                      begin_coarse_chan_idx,
        uintptr_t                      ncoarse_chans
        )
/* Create an array of filenames; free with destroy_filenames()
 */
{
    // Buffer for mwalib error messages
    char error_message[ERROR_MESSAGE_LEN];

    // Check that the number of seconds is non-negative
    if (nseconds <= 0)
    {
        fprintf( stderr, "error: create_filenames: nseconds (= %ld) must be >= 1\n",
                 nseconds );
        exit(EXIT_FAILURE);
    }

    // Check that the number of requested course channels does not exceed the number of available coarse channels
    if (begin_coarse_chan_idx + ncoarse_chans > metafits_metadata->num_metafits_coarse_chans)
    {
        fprintf( stderr, "error: create_filenames: requested coarse channels higher than "
                "what is available\n" );
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the file name list
    char filename[MAX_COMMAND_LENGTH]; // Just the mwalib-generated filename (without the path)
    char **filenames = (char **)malloc( nseconds * ncoarse_chans *sizeof(char *) ); // The full array of filenames, including the paths

    // Allocate memory and write filenames
    unsigned int t_idx, f_idx, second;
    uintptr_t c_idx;
    for (second = 0; second < nseconds; second++)
    {
        t_idx = begin_gps - (metafits_metadata->metafits_timesteps[0].gps_time_ms / 1000) + second;

        for (c_idx = begin_coarse_chan_idx; c_idx < begin_coarse_chan_idx + ncoarse_chans; c_idx++)
        {
            //f_idx = second*ncoarse_chans + c_idx - begin_coarse_chan_idx; // <-- this order seems to break mwalib (v0.9.4)
            f_idx = (c_idx - begin_coarse_chan_idx)*nseconds + second;
            filenames[f_idx] = (char *)malloc( MAX_COMMAND_LENGTH*sizeof(char) );
            if (mwalib_metafits_get_expected_volt_filename(
                        metafits_context,
                        t_idx,
                        c_idx,
                        filename,
                        MAX_COMMAND_LENGTH,
                        error_message,
                        ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
            {
                fprintf( stderr, "error (mwalib): cannot create filenames: %s\n", error_message );
                exit(EXIT_FAILURE);
            }
            sprintf( filenames[f_idx], "%s/%s",
                    datadir, filename );
        }
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


void get_mwalib_metafits_metadata(
        char              *filename,
        MetafitsMetadata **metadata,
        MetafitsContext  **context
        )
{
    char error_message[ERROR_MESSAGE_LEN];
    error_message[0] = '\0'; // <-- Just to avoid a compiler warning about uninitialised variables

    // Create OBS_CONTEXT
    if (mwalib_metafits_context_new2( filename, context, error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( *context, NULL, NULL, metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }
}

void get_mwalib_voltage_metadata(
        VoltageMetadata  **vcs_metadata,
        VoltageContext   **vcs_context,
        MetafitsMetadata **obs_metadata,
        MetafitsContext   *obs_context,
        unsigned long int  begin_gps,
        int                nseconds,
        char              *datadir,
        uintptr_t          coarse_chan_idx,
        int               *ncoarse_chans
        )
/* Create the voltage metadata structs using mwalib's API.
 */
{
    char error_message[ERROR_MESSAGE_LEN];
    error_message[0] = '\0'; // <-- Just to avoid a compiler warning about uninitialised variables

    // If nseconds is an invalid number (<= 0), make it equal to the max possible for this obs
    if (nseconds <= 0)
    {
        nseconds = (*obs_metadata)->num_metafits_timesteps - (begin_gps - (*obs_metadata)->metafits_timesteps[0].gps_time_ms/1000);
    }

    // Ditto for the ncoarse_chans
    if (*ncoarse_chans < 1)
    {
        fprintf( stderr, "error: get_mwalib_voltage_metadata: number of coarse chans must be >1\n" );
        exit(EXIT_FAILURE);
    }

    int nfiles = nseconds * (*ncoarse_chans);

    // Create list of filenames
    char **filenames = create_filenames( obs_context, *obs_metadata, begin_gps, nseconds, datadir, coarse_chan_idx, *ncoarse_chans );

    // Create an mwalib voltage context, voltage metadata, and new obs metadata (now with correct antenna ordering)
    // (MWALIB is expecting a const array, so we will give it one!)
    const char **voltage_files = (const char **)malloc( sizeof(char *) * nfiles );
    int i;
    for (i = 0; i < nfiles; i++)
    {
        voltage_files[i] = filenames[i];
    }

    // Create VCS_CONTEXT
    if (mwalib_voltage_context_new( (*obs_metadata)->metafits_filename, voltage_files, nfiles, vcs_context, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Create VCS_METADATA
    if (mwalib_voltage_metadata_get( *vcs_context, vcs_metadata, error_message, ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot get metadata: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Replace the existing OBS_METADATA with a new one that knows this is a VCS observation
    mwalib_metafits_metadata_free( *obs_metadata );
    if (mwalib_metafits_metadata_get( NULL, NULL, *vcs_context, obs_metadata, error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata from voltage context: %s\n", error_message );
        exit(EXIT_FAILURE);
    }

    // Free memory
    destroy_filenames( filenames, nseconds );
    free( voltage_files );
}


long unsigned int get_relative_gps( MetafitsMetadata *obs_metadata, long int relative_begin )
/* Get the gps second for an observation from a relative value, RELATIVE_BEGIN
 * If RELATIVE_BEGIN >= 0, then return the gps second relative to the "good time"
 * If RELATIVE_BEGIN < 0, then return the gps second relative to the end of the observation
 */
{
    if (relative_begin >= 0)
        return obs_metadata->good_time_gps_ms/1000 + relative_begin;

    // Then relative begin must be negative, so count backwards from end
    return obs_metadata->metafits_timesteps[obs_metadata->num_metafits_timesteps + relative_begin].gps_time_ms/1000;
}

long unsigned int parse_begin_string( MetafitsMetadata *obs_metadata, char *begin_str )
/* Another way to get the beginning gps second
 * If the first character of BEGIN_STR is '+' or '-', then return a relative gps second
 * according to get_relative_gps().
 * Otherwise, parse it as a gps second in its own right.
 */
{
    // If the first character is '+' or '-', treat it as a relative time
    if (begin_str[0] == '+' || begin_str[0] == '-')
        return get_relative_gps( obs_metadata, atol(begin_str) );

    // Otherwise, if the first character is not a number from 0 to 9, generate an error
    if (!(begin_str[0] >= '0' && begin_str[0] <= '9'))
    {
        fprintf( stderr, "error: parse_begin_string: cannot parse '%s' as a valid gps time\n",
                begin_str );
        exit(EXIT_FAILURE);
    }

    // Otherwise, return the number itself
    return atol(begin_str);
}


uintptr_t parse_coarse_chan_string( MetafitsMetadata *obs_metadata, char *begin_coarse_chan_str )
/* Another way to get the coarse channel index (0 to N-1)
 * If the first character of BEGIN_STR is '+' or '-', then return a relative channel idx
 * Otherwise, parse it as a receiver channel number.
 */
{
    if (begin_coarse_chan_str == NULL)
    {
        fprintf( stderr, "error: parse_coarse_chan_string: begin_coarse_chan_str cannot be NULL\n" );
        exit(EXIT_FAILURE);
    }

    // If the first character is '+' or '-', treat it as a relative time
    if (begin_coarse_chan_str[0] == '+')
        return atoi( begin_coarse_chan_str );

    if (begin_coarse_chan_str[0] == '-')
        return obs_metadata->num_metafits_coarse_chans + atoi( begin_coarse_chan_str );

    // Otherwise, if the first character is not a number from 0 to 9, generate an error
    if (!(begin_coarse_chan_str[0] >= '0' && begin_coarse_chan_str[0] <= '9'))
    {
        fprintf( stderr, "error: parse_coarse_chan_string: cannot parse '%s' as a valid coarse channel\n",
                begin_coarse_chan_str );
        exit(EXIT_FAILURE);
    }

    // Otherwise, find the coarse channel with this receiver number
    uintptr_t coarse_chan_idx;
    for (coarse_chan_idx = 0; coarse_chan_idx < obs_metadata->num_metafits_coarse_chans; coarse_chan_idx++)
        if (obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number == (uintptr_t)atoi(begin_coarse_chan_str))
            return coarse_chan_idx;

    // If we got this far, this receiver number doesn't exist in this observation
    fprintf( stderr, "error: parse_coarse_chan_string: receiver coarse chan %s not in this observation\n",
            begin_coarse_chan_str );
    exit(EXIT_FAILURE);

    // Won't ever get here, but return something to avoid compile error/warning
    return 0;
}

