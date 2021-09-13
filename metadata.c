/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <mwalib.h>
#include "metadata.h"

vcsbeam_metadata *init_vcsbeam_metadata(
        char *obs_metafits_filename, char *cal_metafits_filename,
        char *first_coarse_chan_str, int num_coarse_chans_to_process, int coarse_chan_idx_offset,
        char *first_gps_second_str, int num_gps_seconds_to_process, int gps_second_offset,
        char *datadir )
{
    // Allocate memory for the VCSBEAM_METADATA struct
    vcsbeam_metadata *vm = (vcsbeam_metadata *)malloc( sizeof(vcsbeam_metadata) );

    // Get the observation context and metadata
    get_mwalib_metafits_metadata(
            obs_metafits_filename,
            &(vm->obs_metadata),
            &(vm->obs_context) );

    // Convert the (first) chan and gps strings to numbers
    uint32_t first_gps_second = parse_begin_string( vm->obs_metadata, first_gps_second_str ) + gps_second_offset;
    int first_coarse_chan_idx = parse_coarse_chan_string( vm->obs_metadata, first_coarse_chan_str ) + coarse_chan_idx_offset;

    // Get the voltage context and metadata
    vm->num_gps_seconds_to_process  = num_gps_seconds_to_process;
    vm->num_coarse_chans_to_process = num_coarse_chans_to_process;

    get_mwalib_voltage_metadata(
            &(vm->vcs_metadata),
            &(vm->vcs_context),
            &(vm->obs_metadata),
            vm->obs_context,
            first_gps_second,
            num_gps_seconds_to_process,
            datadir,
            first_coarse_chan_idx,
            num_coarse_chans_to_process );

    // Get the calibration context and metadata, if requested
    if (cal_metafits_filename != NULL)
    {
        get_mwalib_metafits_metadata(
                cal_metafits_filename,
                &(vm->cal_metadata),
                &(vm->cal_context) );
    }
    else
    {
        vm->cal_metadata = NULL;
        vm->cal_context = NULL;
    }

    // Construct the gps second and coarse chan idx arrays
    vm->gps_seconds_to_process = (uint32_t *)malloc( num_gps_seconds_to_process * sizeof(uint32_t) );
    vm->coarse_chan_idxs_to_process = (int *)malloc( num_coarse_chans_to_process * sizeof(int) );

    int g;
    for (g = 0; g < num_gps_seconds_to_process; g++)
    {
        vm->gps_seconds_to_process[g] = g + first_gps_second;
    }

    int c;
    for (c = 0; c < num_coarse_chans_to_process; c++)
    {
        vm->coarse_chan_idxs_to_process[c] = c + first_coarse_chan_idx;
    }

    // Set the default output to match the same channelisation as the input
    vm->do_forward_pfb = false;
    vm->do_inverse_pfb = false;

    switch (vm->obs_metadata->mwa_version)
    {
        case VCSLegacyRecombined:
            vm->output_fine_channels = true;
            break;
        case VCSMWAXv2:
            vm->output_coarse_channels = true;
            break;
        default:
            fprintf( stderr, "init_vcsbeam_metadata: error: "
                    "this observation does not appear to be a VCS observation\n" );
            exit(EXIT_FAILURE);
            break;
    }

    // Compute the "shorthand" variables
    vm->sample_rate = vm->vcs_metadata->num_samples_per_voltage_block *
                      vm->vcs_metadata->num_voltage_blocks_per_second;

    vm->bytes_per_second = vm->vcs_metadata->num_voltage_blocks_per_timestep *
                           vm->vcs_metadata->voltage_block_size_bytes;

    // Return the new struct pointer
    return vm;
}


void destroy_vcsbeam_metadata( vcsbeam_metadata *vm )
/* Frees the memory allocated in INIT_VCSBEAM_METADATA
 */
{
    // Free manually created arrays
    free( vm->gps_seconds_to_process );
    free( vm->coarse_chan_idxs_to_process );

    // Free mwalib structs
    mwalib_metafits_metadata_free( vm->obs_metadata );
    mwalib_metafits_context_free( vm->obs_context );
    mwalib_voltage_metadata_free( vm->vcs_metadata );
    mwalib_voltage_context_free( vm->vcs_context );

    if (vm->cal_metadata != NULL)
    {
        mwalib_metafits_metadata_free( vm->cal_metadata );
        mwalib_metafits_context_free( vm->cal_context );
    }

    // Finally, free the struct itself
    free( vm );
}

void set_vcsbeam_fine_output( vcsbeam_metadata *vm, bool switch_on )
/* Turns on/off fine channelised output
 */
{
    vm->output_fine_channels = switch_on;
    if (vm->obs_metadata->mwa_version == VCSMWAXv2)
        vm->do_forward_pfb = switch_on;
}

void set_vcsbeam_coarse_output( vcsbeam_metadata *vm, bool switch_on )
/* Turns on/off coarse channelised output
 */
{
    vm->output_coarse_channels = switch_on;
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
        vm->do_inverse_pfb = switch_on;
}

char **create_filenames(
        const struct MetafitsContext  *metafits_context,
        const struct MetafitsMetadata *metafits_metadata,
        unsigned long int              begin_gps,
        unsigned long int              nseconds,
        char                          *datadir,
        uintptr_t                      begin_coarse_chan_idx,
        uintptr_t                      ncoarse_chans,
        int                           *nfiles
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

    // Work out the number of files
    int timestep_duration_sec      = (metafits_metadata->mwa_version == VCSMWAXv2 ? 8 : 1); // <-- Raise an issue with mwalib (I shouldn't have to do this)
    uint64_t end_gps               = begin_gps + nseconds - 1;
    uint64_t t0_timestep_second    = metafits_metadata->metafits_timesteps[0].gps_time_ms/1000;
    uint64_t start_timestep_second = begin_gps - (begin_gps % timestep_duration_sec);
    uint64_t end_timestep_second   = end_gps - (end_gps % timestep_duration_sec);
    unsigned int ntimesteps        = (end_timestep_second - start_timestep_second)/timestep_duration_sec + 1;
    unsigned int t0_idx            = (start_timestep_second - t0_timestep_second)/timestep_duration_sec;

    *nfiles = ncoarse_chans * ntimesteps;

    // Allocate memory for the file name list
    char filename[MAX_COMMAND_LENGTH]; // Just the mwalib-generated filename (without the path)
    char **filenames = (char **)malloc( *nfiles * ncoarse_chans * sizeof(char *) ); // The full array of filenames, including the paths

    // Allocate memory and write filenames
    unsigned int t_idx, f_idx;
    uintptr_t c_idx;
    for (t_idx = 0; t_idx < ntimesteps; t_idx++)
    {
        for (c_idx = begin_coarse_chan_idx; c_idx < begin_coarse_chan_idx + ncoarse_chans; c_idx++)
        {
            //f_idx = second*ncoarse_chans + c_idx - begin_coarse_chan_idx; // <-- this order seems to break mwalib (v0.9.4)
            f_idx = (c_idx - begin_coarse_chan_idx)*ntimesteps + t_idx;
            filenames[f_idx] = (char *)malloc( MAX_COMMAND_LENGTH );
            memset( filename, 0, MAX_COMMAND_LENGTH );

            if (mwalib_metafits_get_expected_volt_filename(
                        metafits_context,
                        t_idx + t0_idx,
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
        int                ncoarse_chans
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
    if (ncoarse_chans < 1)
    {
        fprintf( stderr, "error: get_mwalib_voltage_metadata: number of coarse chans must be >1\n" );
        exit(EXIT_FAILURE);
    }

    // Create list of filenames
    int nfiles;
    char **filenames = create_filenames( obs_context, *obs_metadata, begin_gps, nseconds, datadir, coarse_chan_idx, ncoarse_chans, &nfiles );

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
    destroy_filenames( filenames, nfiles );
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


int get_num_not_flagged_rf_inputs( vcsbeam_metadata *vm )
{
    int num_not_flagged = 0;
    uintptr_t i;
    for (i = 0; i < vm->obs_metadata->num_rf_inputs; i++)
        if (!(vm->obs_metadata->rf_inputs[i].flagged))
            num_not_flagged++;

    return num_not_flagged;
}
