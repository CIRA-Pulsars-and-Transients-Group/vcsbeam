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

typedef struct vcsbeam_metadata_t
{
    MetafitsContext  *obs_context;    // The mwalib context derived from the target observation's metafits file
    MetafitsMetadata *obs_metadata;   // The mwalib metadata   "      "   "     "         "          "      "

    VoltageContext   *vcs_context;    // The voltage context derived from the available voltage files
    VoltageMetadata  *vcs_metadata;   // The voltage metadata   "      "   "      "        "      "

    MetafitsContext  *cal_context;    // The mwalib context derived from the calibration observation's metafits file
    MetafitsMetadata *cal_metadata;   // The mwalib metadata   "      "   "       "            "         "       "

    int num_coarse_chans_to_process;  // The number of coarse channels to be processed
    int *coarse_chan_idxs_to_process; // A list of the coarse chan idxs to be processed

    int num_gps_seconds_to_process;   // The number of gps seconds to be processed
    uint32_t *gps_seconds_to_process; // A list of the gps seconds to be processed

    bool output_fine_channels;        // Whether to output fine channelised data
    bool output_coarse_channels;      // Whether to output coarse channelised data

    bool do_forward_pfb;              // Whether to perform the forward PFB
    bool do_inverse_pfb;              // Whether to perform the inverse PFB
} vcsbeam_metadata;

/* INIT_VCSBEAM_METADATA
 * =====================
 *
 * Using mwalib, set up the context and metadata structs required to process
 * MWA data. This function sets things up to process a contiguous block of
 * coarse channels and seconds.
 *
 * Inputs:
 *   OBS_METAFITS_FILENAME       - The name of the metafits file for the target observation
 *   CAL_METAFITS_FILENAME       - The name of the metafits file for the associated calibration observation (can be NULL if not required)
 *   FIRST_COARSE_CHAN_STR       - A string representation* of the first coarse channel to be processed (*see below)
 *   NUM_COARSE_CHANS_TO_PROCESS - The number of (contiguous) coarse channels to be processed
 *   COARSE_CHAN_IDX_OFFSET      - Force the processing to begin at a different coarse channel idx
 *   FIRST_GPS_SECOND_STR        - A string representation* of the first gps second to be processed (*see below)
 *   NUM_GPS_SECONDS_TO_PROCESS  - The number of (contiguous) gps seconds to be processed
 *   GPS_SECOND_OFFSET           - Force the processing to begin at a different gps second
 *   DATADIR                     - The folder containing the observation data files
 *
 * Returns:
 *   A pointer to a newly allocated VCSBEAM_METADATA struct
 */
vcsbeam_metadata *init_vcsbeam_metadata(
        char *obs_metafits_filename, char *cal_metafits_filename,
        char *first_coarse_chan_str, int num_coarse_chans_to_process, int coarse_chan_idx_offset,
        char *first_gps_second_str, int num_gps_seconds_to_process, int gps_second_offset,
        char *datadir );


/* DESTROY_VCSBEAM_METADATA
 * ========================
 *
 * Frees the memory allocated in INIT_VCSBEAM_METADATA
 */
void destroy_vcsbeam_metadata( vcsbeam_metadata *vm );

/* SET_VCSBEAM_FINE_OUTPUT & SET_VCSBEAM_COARSE_OUTPUT
 * =======================   =========================
 *
 * Turns on/off fine/coarse channelised output
 */
void set_vcsbeam_fine_output( vcsbeam_metadata *vm, bool switch_on );
void set_vcsbeam_coarse_output( vcsbeam_metadata *vm, bool switch_on );

// OTHER AUXILIARY FUNCTIONS

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
        int                ncoarse_chans
        );

long unsigned int get_relative_gps( MetafitsMetadata *obs_metadata, long int relative_begin );
long unsigned int parse_begin_string( MetafitsMetadata *obs_metadata, char *begin_str );
uintptr_t parse_coarse_chan_string( MetafitsMetadata *obs_metadata, char *begin_coarse_chan_str );

int get_num_not_flagged_rf_inputs( vcsbeam_metadata *vm );
#endif
