/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __VCSBEAM_H__
#define __VCSBEAM_H__

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>

#include <mpi.h>

#include <mwalib.h>
#include <star/pal.h>
#include <star/palmac.h>
#include <psrfits.h>


#define VCSBEAM_VERSION  "v2.14.3_8d43a51"
#define RUNTIME_DIR      "/home/smcsweeney/bin/vcsbeam_runtime"
#define HYPERBEAM_HDF5   "/opt/mwa_hyperbeam/mwa_full_embedded_element_pattern.h5"


/*******************
 *                 *
 *     MWALIB      *
 *                 *
 *******************/

#define ERROR_MESSAGE_LEN  1024
#define MAX_COMMAND_LENGTH 1024

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



/*******************
 *                 *
 *   PERFORMANCE   *
 *                 *
 *******************/

#define PERFORMANCE_MAX_NUM_STOPWATCHES  16
#define PERFORMANCE_MAX_START_STOP       4096

#define PERFORMANCE_NO_STOPWATCH_FOUND  -1

#define PERFORMANCE_NO_MPI  -1

typedef struct logger_stopwatch_t
{
    char    *name;
    char    *description;
    double   values[PERFORMANCE_MAX_START_STOP];
    int      nstart_stops;
    bool     running;
    double   total, total_sq; // For calculation of stats
} logger_stopwatch;

typedef struct logger_t
{
    double            begintime;
    FILE             *fout;
    logger_stopwatch  stopwatches[PERFORMANCE_MAX_NUM_STOPWATCHES];
    int               nstopwatches;
    int               world_rank;
} logger;

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************
 * Functions for memory management and initialisation *
 ******************************************************/

logger *create_logger( FILE * fout, int world_rank );
void destroy_logger( logger *log );

/*******************************************************
 * Functions for manipulating user-defined stopwatches *
 *******************************************************/

void logger_add_stopwatch( logger *log, const char *stopwatch_name, const char *description );
void logger_start_stopwatch( logger *log, const char *stopwatch_name, bool print_description );
void logger_stop_stopwatch( logger *log, const char *stopwatch_name );

/********************************************************************
 * Functions for writing out different kinds of messages to the log *
 ********************************************************************/

void logger_timed_message( logger *log, const char *message );
void logger_message( logger *log, const char *message );
void logger_stopwatch_report_stats( logger *log, const char *stopwatch_name );
void logger_report_all_stats( logger *log );

#ifdef __cplusplus
} // End extern "C"
#endif

/*******************
 *                 *
 *     FILTER      *
 *                 *
 *******************/

typedef enum filter_type_t
{
    ANALYSIS_FILTER,
    SYNTHESIS_FILTER
} filter_type;

typedef struct pfb_filter_t
{
    double          *coeffs;
    int              ncoeffs;
    int              ntaps;
    int              nchans; // = size/ntaps
    cuDoubleComplex *twiddles; // twiddle factors
    filter_type      type;
} pfb_filter;



#ifdef __cplusplus
extern "C" {
#endif

/* ROOTS_OF_UNITY
 * ==============
 *
 * Creates a complex-valued array containing the N roots of unity.
 * The caller should free this memory (via free()).
 */
cuDoubleComplex *roots_of_unity( int N );



/* LOAD_FILTER_COEFFICIENTS
 * ========================
 *
 * Load a set of filter coefficients
 * Inputs:
 *   FILTERNAME - string specifying a filter. There should be a corresponding
 *                file in the RUNTIME_DIR called FILTERNAME.dat.
 *   TYPE       - Whether it's an ANALYSIS_FILTER or a SYNTHESIS_FILTER
 *   NCHANS     - the number of channels that this filter will be applied to.
 *                For both ANALYSIS and SYNTHESIS filters, this should be
 *                the number of ANALYSIS channels.
 * Outputs:
 *   [return value] - pointer to the newly allocated struct containing the
 *                    coefficients. This should be freed using
 *                    free_pfb_filter()
 */
pfb_filter *load_filter_coefficients( char *filtername, filter_type type, int nchans );



/* FREE_PFB_FILTER
 * ===============
 *
 * Free the memory allocated in load_filter_coefficients()
 */
void free_pfb_filter( pfb_filter *filter );


#ifdef __cplusplus
} // End extern "C"
#endif


/*******************
 *                 *
 *    METADATA     *
 *                 *
 *******************/


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

    // Some "shorthand" variables
    // These can be worked out from the other fields, but are computed here
    // for convenience.
    int sample_rate;                  // Number of samples per second
    int bytes_per_second;             // Bytes per second of data
} vcsbeam_metadata;


#ifdef __cplusplus
extern "C" {
#endif

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
        uintptr_t                      ncoarse_chans,
        int                           *nfiles
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

#ifdef __cplusplus
} // End extern "C"
#endif

/*******************
 *                 *
 *    GEOMETRY     *
 *                 *
 *******************/

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

struct beam_geom {
    double mean_ra;
    double mean_dec;
    double az;
    double el;
    double lmst;
    double fracmjd;
    double intmjd;
    double unit_N;
    double unit_E;
    double unit_H;
};


/********************
 *                  *
 *   BEAM_PSRFITS   *
 *                  *
 ********************/

typedef struct mpi_psrfits_t
{
    MPI_Datatype    coarse_chan_spectrum;
    MPI_Datatype    total_spectrum_type;
    MPI_Datatype    spliced_type;

    MPI_Request     request_data;
    MPI_Request     request_offsets;
    MPI_Request     request_scales;

    int             ncoarse_chans;

    struct psrfits  coarse_chan_pf;
    struct psrfits  spliced_pf;

    int             writer_id;
} mpi_psrfits;



/********************
 *                  *
 *   FORWARD PFB    *
 *                  *
 ********************/

/* The final step in the forward PFB version that emulates the FPGA
 * implementation packs the data into (4+4)-bit complex samples. The
 * following macros collectively achieve this.
 */
#define CLIP(x,max)        ((x) < -(max)   ? -(max)   : \
                            (x) >  (max)-1 ?  (max)-1 : (x))
#define INT_TO_UINT8(x)    ((x) < 0 ? CLIP(x,8) + 0x10 : CLIP(x,8))
#define DEMOTE(x)  (INT_TO_UINT8((int)round(x)))
#define PACK_NIBBLES(r,i)  ((DEMOTE(i) << 4) + DEMOTE(r))

typedef enum pfb_result_t
{
    PFB_SUCCESS,
    PFB_END_OF_GPSTIMES
} pfb_result;

typedef enum pfb_flags_t
{
    PFB_MALLOC_HOST_INPUT    = 0x01,
    PFB_MALLOC_HOST_OUTPUT   = 0x02,
    PFB_MALLOC_DEVICE_INPUT  = 0x04,
    PFB_MALLOC_DEVICE_OUTPUT = 0x08,
    PFB_MALLOC_ALL           = 0x0F
} pfb_flags;

typedef struct forward_pfb_t
{
    vcsbeam_metadata *vm;                     // All the necessary metadata

    char2            *htr_data;               // The input data, as obtained via mwalib from coarse channelised data
    char2            *d_htr_data;             // Same as above, on device

    uint8_t          *vcs_data;               // The output data, fine channelised and packed into the VCS recombined format
    uint8_t          *d_vcs_data;             // Same as above, on device

    size_t            htr_size;               // The size (in bytes) of htr_data
    size_t            vcs_size;               // The size (in bytes) of vcs_data

    size_t            char2s_per_second;      // The number of char2's in one second of HTR data
    size_t            bytes_per_block;        // The number of bytes in one "voltage block" of HTR data

    size_t            current_gps_idx;        // Which gps second (in vm) to read next

    int               ninputs_per_cufft_batch; // Necessary because one can't do 2560000 batches, apparently
    int               cufft_batch_size;

    cuFloatComplex   *d_weighted_overlap_add; // A "temporary" array on the device for mid-calculation product
    size_t            weighted_overlap_add_size; // The size (in bytes) of d_weighted_overlap_add

    int              *filter_coeffs;          // The filter to be applied **WARNING! Filter will be typecast to int!!**
    int              *d_filter_coeffs;        // As above, on the device

    int               nspectra;               // The number of spectra to generate
    int               M;                      // The "stride" of the PFB (setting M=K means critically sampled)
    int               K;                      // The number of channels
    int               I;                      // The number of RF inputs
    int               P;                      // The number of taps

    cufftHandle       plan;                   // The cuFFT plan for performing the FFT part of the forward PFB

    bool              read_locked;            // For multi-threading: lock for reading
} forward_pfb;


#ifdef __cplusplus
extern "C" {
#endif

forward_pfb *init_forward_pfb( vcsbeam_metadata *vm, pfb_filter *filter, int M, pfb_flags flags );

void free_forward_pfb( forward_pfb *fpfb );

pfb_result forward_pfb_read_next_second( forward_pfb *fpfb );
    
void cu_forward_pfb_fpga_version( forward_pfb *fpfb, bool copy_result_to_host, logger *log );

#ifdef __cplusplus
}
#endif




/********************
 *                  *
 *   BACKWARD PFB   *
 *                  *
 ********************/

struct gpu_ipfb_arrays
{
    int ntaps;
    int in_size;
    int ft_size;
    int out_size;
    float *in_real,   *in_imag;
    float *ft_real,   *ft_imag;
    float *d_in_real, *d_in_imag;
    float *d_ft_real, *d_ft_imag;
    float *d_out;
};


#ifdef __cplusplus
extern "C" {
#endif

void cu_invert_pfb( cuDoubleComplex ****detected_beam, int file_no,
                        int npointing, int nsamples, int nchan, int npol, int sizeof_buffer,
                        struct gpu_ipfb_arrays *g, float *data_buffer_uvdif );

void cu_load_ipfb_filter( pfb_filter *filter, struct gpu_ipfb_arrays *g );

void malloc_ipfb( struct gpu_ipfb_arrays *g, pfb_filter *filter, int nsamples,
        int npol, int npointing );

void free_ipfb( struct gpu_ipfb_arrays *g );


#ifdef __cplusplus
}
#endif




/********************
 *                  *
 *    FORM BEAM     *
 *                  *
 ********************/

/* Converting from 4+4 complex to full-blown complex doubles */

#define REAL_NIBBLE_TO_UINT8(X)  ((X) & 0xf)
#define IMAG_NIBBLE_TO_UINT8(X)  (((X) >> 4) & 0xf)
#define UINT8_TO_INT(X)          ((X) >= 0x8 ? (signed int)(X) - 0x10 : (signed int)(X))
#define RE_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))))
#define IM_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X))))
#define UCMPLX4_TO_CMPLX_FLT(X)  (make_cuDoubleComplex((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))), \
                                         (float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X)))))
#define DETECT(X)                (cuCreal(cuCmul(X,cuConj(X))))


/* structure for managing data arrays to be allocated on both host and device */
struct gpu_formbeam_arrays
{
    size_t coh_size;
    size_t data_size;
    size_t Bd_size;
    size_t J_size;
    size_t JD_size;
    size_t pol_idxs_size;
    cuDoubleComplex *J, *d_J;
    cuDoubleComplex *Bd, *d_Bd;
    cuDoubleComplex *JDx, *d_JDx;
    cuDoubleComplex *JDy, *d_JDy;
    uint8_t *d_data;
    float   *d_coh;
    uint32_t *polX_idxs, *d_polX_idxs;
    uint32_t *polY_idxs, *d_polY_idxs;
};

#ifdef __cplusplus
extern "C" {
#endif

void malloc_formbeam( struct gpu_formbeam_arrays *g, vcsbeam_metadata *vm,
                      int *nchunk, float gpu_mem_gb, int outpol_coh,
                      int npointing, logger *log );
void free_formbeam( struct gpu_formbeam_arrays *g );



void cu_form_incoh_beam(
        uint8_t *data, uint8_t *d_data, size_t data_size,
        float *d_incoh,
        unsigned int nsample, int nchan, int ninput,
        float *offsets, float *d_offsets,
        float *scales, float *d_scales,
        uint8_t *Iscaled, uint8_t *d_Iscaled, size_t Iscaled_size
        );


void cu_form_beam( uint8_t *data, unsigned int sample_rate, cuDoubleComplex *d_phi,
                   int file_no,
                   int npointing, int nstation, int nchan,
                   int npol, double invw, struct gpu_formbeam_arrays *g,
                   cuDoubleComplex ****detected_beam, float *coh,
                   cudaStream_t *streams, int nchunk,
                   mpi_psrfits *mpfs );

float *create_pinned_data_buffer_psrfits( size_t size );

float *create_pinned_data_buffer_vdif( size_t size );

void cu_upload_pol_idx_lists( struct gpu_formbeam_arrays *g );

cuDoubleComplex ****create_invJi( int nstation, int nchan, int pol );
void              destroy_invJi( cuDoubleComplex ****array, int nstation, int nchan, int npol );

cuDoubleComplex ****create_detected_beam( int npointing, int nsamples, int nchan, int npol );
void              destroy_detected_beam( cuDoubleComplex ****array, int npointing,
                                         int nsamples, int nchan );

void allocate_input_output_arrays( void **data, void **d_data, size_t size );
void free_input_output_arrays( void *data, void *d_data );

void flatten_bandpass( int nstep, int nchan, int npol, void *data);

#ifdef __cplusplus
}
#endif

#endif
