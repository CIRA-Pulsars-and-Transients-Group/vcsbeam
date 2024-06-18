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
#include <inttypes.h>

#include "gpu_includes.h"
#include "gpu_macros.h"
#include "gpu_fft.hpp"
#include <mpi.h>

#include <mwalib.h>
#include <star/pal.h>
#include <star/palmac.h>
#include <psrfits.h>
#include <vdifio.h>
#include <mwa_hyperbeam.h>



#define VCSBEAM_VERSION  "v4.2.95_26a66c7"
#define RUNTIME_DIR      "/usr/local/bin/vcsbeam_runtime"
/* #undef HYPERBEAM_HDF5 */


/* Boilerplate CUDA code for error checking */

/*#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)*/


/*******************
 *                 *
 *     MWALIB      *
 *                 *
 *******************/

#define ERROR_MESSAGE_LEN  1024
#define MAX_COMMAND_LENGTH 1024

/* The following is a patch to provide constants which are not available in mwalib */

#ifndef SPEED_OF_LIGHT_IN_VACUUM_M_PER_S
#define SPEED_OF_LIGHT_IN_VACUUM_M_PER_S MWALIB_SPEED_OF_LIGHT_IN_VACUUM_M_PER_S
#endif

#ifndef MWA_LATITUDE_RADIANS
#define MWA_LATITUDE_RADIANS MWALIB_MWA_LATITUDE_RADIANS
#endif

#ifndef MWA_LONGITUDE_RADIANS
#define MWA_LONGITUDE_RADIANS MWALIB_MWA_LONGITUDE_RADIANS
#endif

#ifndef MWA_ALTITUDE_METRES
#define MWA_ALTITUDE_METRES MWALIB_MWA_ALTITUDE_METRES
#endif

/* End replacement constants */

/* Get angle constants from PAL, if available;
 * otherwise construct them from math.h
 */
#ifdef PAL_FOUND
#define PIBY2  PAL__DPIBY2
#define H2R    PAL__DH2R
#define R2H    PAL__DR2H
#define D2R    PAL__DD2R
#define R2D    PAL__DR2D
#else
#define PIBY2  (0.5*M_PI)
#define H2R    (M_PI/12.0)
#define R2H    (12.0/M_PI)
#define D2R    (M_PI/180.0)
#define R2D    (180.0/M_PI)
#endif


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
 *   CALIBRATION   *
 *                 *
 *******************/

// Calibration solution types
#define CAL_NONE      0
#define CAL_RTS       1
#define CAL_OFFRINGA  2

#define OFFRINGA_HEADER_SIZE_BYTES 48
#define JONES_SIZE_BYTES           64  /* (4 elements) * (2 re/im) * (8 bytes per double) */

#define NDBL_PER_JONES  8
#define BANDPASS_ROWS_PER_ANT  8

#define CAL_BUFSIZE    4096


typedef struct calibration_t
{
    char  *metafits;           // Filename of the metafits file
    char  *caldir;             // Location of calibration data
    int    cal_type;           // Either RTS or OFFRINGA
    char  *ref_ant;            // Reference antenna for calibration phases
    bool   keep_cross_terms;   // Include PQ and QP of calibration Jones matrices
    double phase_offset;       // Rotate the phase of Y by m*freq + c, where
    double phase_slope;        //   m = phase_slope (rad/Hz)
                               //   c = phase_offset (rad)
    bool   use_bandpass;       // Use the Bandpass solutions
    bool   custom_pq_correction; // Set to true if phase_offset and phase_slope are to be applied
    char  *flags_file;         // Name of file containing custom flagged tiles
    int    nflags;             // The size of flagged_tilenames
    char **flagged_tilenames;  // a list of tilenames to be flagged
} calibration;


#ifdef __cplusplus
extern "C" {
#endif


void read_dijones_file( gpuDoubleComplex **Dd, gpuDoubleComplex *A, double *amp, uintptr_t nant, char *fname );
void read_bandpass_file( gpuDoubleComplex ***Jm, gpuDoubleComplex ***Jf,
        MetafitsMetadata *cal_metadata, char *filename );

void remove_reference_phase( gpuDoubleComplex *J, gpuDoubleComplex *Jref );
void zero_PQ_and_QP( gpuDoubleComplex *J );

void parse_calibration_correction_file( uint32_t gpstime, calibration *cal );

bool tilename_is_flagged( char *tilename, calibration *cal );

void init_calibration( calibration *cal );
void free_calibration( calibration *cal );

#ifdef __cplusplus
}
#endif


/******************
 *                *
 *     JONES      *
 *                *
 ******************/

#define NSTOKES  4

/* Calculating array indices for various flavours of Jones matrix arrays */

#define v_IDX(s,c,i,nc,ni)    ((s) * (ni)*(nc) + \
                               (c) * (ni)      + \
                               (i))

#define VSPVB  64000  /* "Vcsmwax Samples Per Voltage Block"
                         (i.e. the size of the TIME dimension in a voltage block */
#define vMWAX_IDX(s,i,ni) (((s/VSPVB)*(ni) + (i))*VSPVB + (s%VSPVB))

#define D_IDX(a,c,p1,p2,nc,npol)   ((a)  * (npol)*(npol)*(nc)        + \
                                    (c)  * (npol)*(npol)             + \
                                    (p1) * (npol)                    + \
                                    (p2))

#define J_IDX(p,a,c,p1,p2,nant,nc,npol)   ((p)  * (nant)*(npol)*(npol)*(nc) + \
                                           (a)  * (npol)*(npol)*(nc)        + \
                                           (c)  * (npol)*(npol)             + \
                                           (p1) * (npol)                    + \
                                           (p2))

#define Jv_IDX(p,s,c,a,ns,nc,na)  ((p) * (na)*(nc)*(ns) + \
                                   (s) * (na)*(nc)      + \
                                   (c) * (na)           + \
                                   (a))

#define B_IDX(p,s,c,pol,ns,nc,npol) ((p)  * (npol)*(nc)*(ns)   + \
                                     (s)  * (npol)*(nc)        + \
                                     (c)  * (npol)             + \
                                     (pol))

#define C_IDX(p,s,st,c,ns,nst,nc)  ((p)  * ((nc)*(nst)*(ns)) + \
                                    (s)  * ((nc)*(nst))      + \
                                    (st) *  (nc)               + \
                                    (c))

#define I_IDX(s,c,nc)          ((s)*(nc) + (c))

#ifdef __cplusplus
extern "C" {
#endif

/**** MATRIX OPERATIONS ****/

void cp2x2( gpuDoubleComplex *Min, gpuDoubleComplex *Mout );
void inv2x2( gpuDoubleComplex *Min, gpuDoubleComplex *Mout );
void inv2x2d( double *Min, double *Mout );
void inv2x2S( gpuDoubleComplex *Min, gpuDoubleComplex *Mout );
void mult2x2d( gpuDoubleComplex *M1, gpuDoubleComplex *M2, gpuDoubleComplex *Mout );
void mult2x2d_RxC( double *M1, gpuDoubleComplex *M2, gpuDoubleComplex *Mout );
void mult2x2d_CxR( gpuDoubleComplex *M1, double *M2, gpuDoubleComplex *Mout );
void conj2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout );
double norm2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout );
void reverse2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout );
void swaprows2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout );
void swapcols2x2( gpuDoubleComplex *M, gpuDoubleComplex *Mout );
bool is2x2zero( gpuDoubleComplex *M );
void calc_hermitian( gpuDoubleComplex *M, gpuDoubleComplex *H );
void calc_coherency_matrix( gpuDoubleComplex *M, gpuDoubleComplex *c );
void fprintf_complex_matrix( FILE *fout, gpuDoubleComplex *M );

#ifdef __cplusplus
}
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
    int              nchans; // = ncoeffs/ntaps
    gpuDoubleComplex *twiddles; // twiddle factors
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
gpuDoubleComplex *roots_of_unity( int N );



#ifdef __cplusplus
} // End extern "C"
#endif


/*******************
 *                 *
 *    METADATA     *
 *                 *
 *******************/

typedef enum vcsbeam_datatype_t
{
    VM_INT4,
    VM_DBL
} vcsbeam_datatype;

typedef enum vm_error_t
{
    VM_SUCCESS,
    VM_END_OF_DATA,
    VM_READ_BUFFER_NOT_SET,
    VM_READ_BUFFER_LOCKED
} vm_error;

typedef struct geometric_delays_t {
    gpuDoubleComplex   *phi;
    gpuDoubleComplex   *d_phi;
    uintptr_t          npointings;
    uintptr_t          nant;
    uintptr_t          nchan;
    double            *chan_freqs_hz;
    MetafitsMetadata  *obs_metadata;
} geometric_delays;

typedef struct primary_beam_t
{
    gpuDoubleComplex  *B;
    FEEBeam          *beam;
    uint32_t        **delays;
    double          **amps;
    uintptr_t         npointings;
    uintptr_t         nant;
    uintptr_t         npol;
    uint32_t          freq_hz;
    MetafitsMetadata *obs_metadata;
} primary_beam;

typedef struct beam_geom_t {
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
} beam_geom;


typedef struct mpi_psrfits_t
{
    MPI_Datatype    coarse_chan_spectrum;
    MPI_Datatype    total_spectrum_type;
    MPI_Datatype    spliced_type;

    MPI_Datatype    coarse_chan_OFFS;
    MPI_Datatype    total_OFFS;
    MPI_Datatype    spliced_OFFS;

    MPI_Datatype    coarse_chan_SCL;
    MPI_Datatype    total_SCL;
    MPI_Datatype    spliced_SCL;
    
    MPI_Request     request_data;
    MPI_Request     request_offsets;
    MPI_Request     request_scales;

    struct psrfits  coarse_chan_pf;
    struct psrfits  spliced_pf;

    int             writer_id;
} mpi_psrfits;

typedef struct host_buffer_t
{
    void   *buffer;        // The beginning of the whole buffer
    size_t  buffer_size;   // The size (in bytes) of the whole buffer

    void   *read_ptr;      // The pointer of where to read into
    size_t  read_size;     // The number of bytes to read in

    void   *copy_from_ptr; // The pointer of where to copy from
    void   *copy_to_ptr;   // The pointer of where to copy to
    size_t  copy_size;     // The number of bytes to copy

    bool    locked;        // Whether reading/copying is allowed
} host_buffer;


typedef struct device_buffer_t
{
    void   *buffer;
    size_t  buffer_size;
} device_buffer;

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
#define PACK_NIBBLES(r,i)  ((DEMOTE(r) << 4) + DEMOTE(i))

typedef enum pfb_flags_t
{
    PFB_MALLOC_HOST_INPUT    = 0x01,
    PFB_MALLOC_HOST_OUTPUT   = 0x02,
    PFB_MALLOC_DEVICE_INPUT  = 0x04,
    PFB_MALLOC_DEVICE_OUTPUT = 0x08,
    PFB_MALLOC_ALL           = 0x0F,

    PFB_TYPE_MASK            = 0xF0, // Next four bytes used for listing different output data types
    PFB_COMPLEX_INT4         = 0x10,
    PFB_COMPLEX_FLOAT64      = 0x20,

    PFB_IMAG_PART_FIRST      = 0x100,

    PFB_EMULATE_FPGA         = 0x200,

    // Some summary flag settings:
    PFB_SMART                = 0x21D, // = PFB_MALLOC_HOST_INPUT | PFB_MALLOC_DEVICE_INPUT | PFB_MALLOC_DEVICE_OUTPUT |
                                      //   PFB_EMULATE_FPGA | PFB_COMPLEX_INT4
    PFB_FULL_PRECISION       = 0x2D   // = PFB_MALLOC_HOST_INPUT | PFB_MALLOC_DEVICE_INPUT | PFB_MALLOC_DEVICE_OUTPUT |
                                      //   PFB_COMPLEX_FLOAT64
} pfb_flags;

typedef struct forward_pfb_t
{
    char2            *d_htr_data;             // Same as above, on device

    void             *vcs_data;               // The output data, fine channelised and packed into the VCS recombined format
    void             *d_vcs_data;             // Same as above, on device

    size_t            d_htr_size;             // The size (in bytes) of d_htr_data
    size_t            htr_stride;             // Stride for chunks
    size_t            vcs_size;               // The size (in bytes) of vcs_data
    size_t            d_vcs_size;             // The size (in bytes) of d_vcs_data
    size_t            vcs_stride;             // Stride for chunks

    size_t            char2s_per_second;      // The number of char2's in one second of HTR data
    size_t            bytes_per_block;        // The number of bytes in one "voltage block" of HTR data
    size_t            chunks_per_second;      // Split each second into this many chunks on the GPU
    size_t            chunk;                  // The current chunk being/about to be processed

    int               ninputs_per_cufft_batch; // Necessary because one can't do 2560000 batches, apparently
    int               cufft_batch_size;

    gpuFloatComplex   *d_weighted_overlap_add; // A "temporary" array on the device for mid-calculation product
    size_t            weighted_overlap_add_size; // The size (in bytes) of d_weighted_overlap_add

    int              *filter_coeffs;          // The filter to be applied **WARNING! Filter will be typecast to int!!**
    int              *d_filter_coeffs;        // As above, on the device

    int               nspectra;               // The number of spectra to generate
    int               nspectra_per_chunk;     // The number of spectra per chunk
    int               M;                      // The "stride" of the PFB (setting M=K means critically sampled)
    int               K;                      // The number of channels
    int               I;                      // The number of RF inputs
    int               P;                      // The number of taps

    int              *i_output_idx;           // The idxs for where to put each RF input in the output
    int              *d_i_output_idx;         // (The MWAX files are ordered by "Antenna", whereas the legacy
                                              // VCS files are ordered in their own special order)

    pfb_flags         flags;                  // See pfb_flags enum above for options
    gpufftHandle       plan;                   // The cuFFT plan for performing the FFT part of the forward PFB
} forward_pfb;


typedef struct vcsbeam_context_t
{
    // MPI
    bool use_mpi;                     // Toggle MPI usage
    int mpi_size;                     // MPI world size
    int mpi_rank;                     // MPI process rank
    int ncoarse_chans;
    int coarse_chan_idx;
    int writer;                       // The rank of the process responsible for writing output files

    // Observation metadata
    MetafitsContext  *obs_context;    // The mwalib context derived from the target observation's metafits file
    MetafitsMetadata *obs_metadata;   // The mwalib metadata   "      "   "     "         "          "      "

    VoltageContext   *vcs_context;    // The voltage context derived from the available voltage files
    VoltageMetadata  *vcs_metadata;   // The voltage metadata   "      "   "      "        "      "

    MetafitsContext  *cal_context;    // The mwalib context derived from the calibration observation's metafits file
    MetafitsMetadata *cal_metadata;   // The mwalib metadata   "      "   "       "            "         "       "

    MetafitsContext  *obs_context_legacy; // The same as obs_context, but forced to be legacy
    MetafitsMetadata *obs_metadata_legacy; // The same as obs_context, but forced to be legacy

    char *datadir;                    // The directory containing the input data
    char **filenames;                 // The input filenames
    int nfiles;                       // The number of filenames
    unsigned int nfiletimes;          // The number of "file" timesteps
    unsigned int seconds_per_file;    // The number of seconds in one file

    // Calibration
    calibration cal;                  // Calibration specification

    // Tied-array beam pointings
    double *ras_hours;                // Array of pointing RAs (in decimal hours)
    double *decs_degs;                // Array of pointing degs (in decimal degrees)

    // PFB
    pfb_filter *analysis_filter;      // Filters for PFB analysis
    pfb_filter *synth_filter;         // Filter for PFB synthesis

    forward_pfb *fpfb;                // Forward PFB

    // Data buffers
    host_buffer *v;                   // The buffer for the input data on host
    void *d_v;                        // The buffer for the input data on device
    uintptr_t v_size_bytes;           // The size of data in bytes (currently always = bytes_per_second)
    uintptr_t d_v_size_bytes;         // The size of d_data in bytes (depends on number of "chunks")

    void *S, *d_S;                    // The buffers for the detected (full Stokes) coherent beam on host/device
    uintptr_t S_size_bytes;           // The size of S in bytes
    uintptr_t d_S_size_bytes;         // The size of d_S in bytes

    gpuDoubleComplex *e, *d_e;         // The buffers for the beamformed voltages on host/device
    uintptr_t e_size_bytes;           // The size of S in bytes
    uintptr_t d_e_size_bytes;         // The size of d_S in bytes

    gpuDoubleComplex *J, *d_J;         // The buffers for Jones matrices on host/device
    uintptr_t J_size_bytes;           // The size of J in bytes
    uintptr_t d_J_size_bytes;         // The size of d_J in bytes

    gpuDoubleComplex *Jv_P, *d_Jv_P;   // The buffers for Jones-corrected voltages on host/device
    gpuDoubleComplex *Jv_Q, *d_Jv_Q;   // The buffers for Jones-corrected voltages on host/device
    uintptr_t Jv_size_bytes;          // The size of Jv in bytes
    uintptr_t d_Jv_size_bytes;        // The size of d_Jv in bytes

    gpuDoubleComplex *D, *d_D;         // The buffers for calibration solutions on host/device
    uintptr_t D_size_bytes;           // The size of D in bytes
    uintptr_t d_D_size_bytes;         // The size of d_D in bytes

    uint32_t *polP_idxs, *d_polP_idxs; // List of indices for VCS-ordered data
    uint32_t *polQ_idxs, *d_polQ_idxs;
    uintptr_t pol_idxs_size_bytes;     // The size of (each of) the P/Q idxs (host) arrays
    uintptr_t d_pol_idxs_size_bytes;   // The size of (each of) the P/Q idxs (device) arrays

    geometric_delays gdelays;         // Geometric delays due to tile layout
    primary_beam pb;                  // Primary beam calculations

    unsigned int npointing;           // Number of requested tied array beam pointings

    uintptr_t max_gpu_mem_bytes;      // The maximum allowed GPU memory to use (in bytes)
    uint32_t chunks_per_second;       // The number of chunks to process on device per second of data
                                      // (data_size_bytes = d_data_size_bytes * chunks_per_second)

    int num_coarse_chans_to_process;  // The number of coarse channels to be processed
    int *coarse_chan_idxs_to_process; // A list of the coarse chan idxs to be processed
    int *cal_coarse_chan_idxs_to_process; // A list of the calibration coarse chan idxs to be processed

    int num_gps_seconds_to_process;   // The number of gps seconds to be processed
    uint32_t *gps_seconds_to_process; // A list of the gps seconds to be processed

    bool output_fine_channels;        // Whether to output fine channelised data
    bool output_coarse_channels;      // Whether to output coarse channelised data

    size_t current_gps_idx;           // Which gps second to read next

    bool do_forward_pfb;              // Whether to perform the forward PFB
    bool do_inverse_pfb;              // Whether to perform the inverse PFB

    // Some "shorthand" variables
    // These can be worked out from the other fields, but are computed here
    // for convenience.
    int sample_rate;                  // Number of samples per second
    int fine_sample_rate;             // Number of samples per second for fine channelised data
    int bytes_per_second;             // Bytes per second of data
    int num_not_flagged;              // Number of RF inputs that are not flagged
    int nchan;                        // Number of channels
    int nfine_chan;                   // Number of fine channels

    // Output number of stokes parameters needed
    int out_nstokes;

    // VDIF output
    vdif_header      vhdr;
    struct vdifinfo *vf;

    // Beam output statistics
    float *d_offsets, *offsets;
    float *d_scales, *scales;
    uint8_t *d_Cscaled, *Cscaled;

    size_t offsets_size;
    size_t scales_size;
    size_t Cscaled_size;

    // *** FOR INTERNAL USE ONLY ***
    int chunk_to_load;                // Which chunk number to load to gpu next
    vcsbeam_datatype datatype;        // What format the input data samples are in

    gpuStream_t *streams;            // CUDA streams used for forming multiple tied-array beams

    logger *log;                      // Used for log messages
    char log_message[MAX_COMMAND_LENGTH];

    char error_message[ERROR_MESSAGE_LEN]; // MWALIB error message buffer
} vcsbeam_context;

#define NO_ANTENNA_FOUND  -1

#ifdef __cplusplus
extern "C" {
#endif

void vmCheckError( vm_error err );

vcsbeam_context *vmInit( bool use_mpi );

/**    
 * vmBindObsData
 * =============
 *
 * Using mwalib, set up the context and metadata structs required to process
 * MWA data. This function sets things up to process a contiguous block of
 * coarse channels and seconds.
 *
 * Inputs:
 *   VM                          - VCSBeam Context
 *   FIRST_COARSE_CHAN_STR       - A string representation* of the first coarse channel to be processed (*see below)
 *   NUM_COARSE_CHANS_TO_PROCESS - The number of (contiguous) coarse channels to be processed
 *   COARSE_CHAN_IDX_OFFSET      - Force the processing to begin at a different coarse channel idx
 *   FIRST_GPS_SECOND_STR        - A string representation* of the first gps second to be processed (*see below)
 *   NUM_GPS_SECONDS_TO_PROCESS  - The number of (contiguous) gps seconds to be processed
 *   GPS_SECOND_OFFSET           - Force the processing to begin at a different gps second
 *   DATADIR                     - The folder containing the observation data files
 */
void vmBindObsData(
        vcsbeam_context *vm,
        char *first_coarse_chan_str, int num_coarse_chans_to_process, int coarse_chan_idx_offset,
        char *first_gps_second_str, int num_gps_seconds_to_process, int gps_second_offset,
        char *datadir );

void vmBindCalibrationData( vcsbeam_context *vm,
        char   *caldir,
        int     cal_type,
        bool    use_bandpass,  // only relevant for RTS observations
        char   *flags_file );

void vmReadCalibration( vcsbeam_context *vm );

void vmPrintTitle( vcsbeam_context *vm, const char *title );

/* DESTROY_VCSBEAM_METADATA
 * ========================
 *
 * Frees the memory allocated in INIT_VCSBEAM_METADATA
 */
void destroy_vcsbeam_context( vcsbeam_context *vm );

/**
 * vmSetOutputChannelisation
 * =========================
 *
 * Turns on/off fine/coarse channelised output, and deduces whether the
 * forward/inverse pfb will be needed, based on the input channelisation
 */
void vmSetOutputChannelisation( vcsbeam_context *vm, bool out_fine, bool out_coarse );

// OTHER AUXILIARY FUNCTIONS

void vmMallocVHost( vcsbeam_context *vm );
void vmMallocVDevice( vcsbeam_context *vm );
void vmFreeVHost( vcsbeam_context *vm );
void vmFreeVDevice( vcsbeam_context *vm );

void vmMallocJVHost( vcsbeam_context *vm );
void vmMallocJVDevice( vcsbeam_context *vm );
void vmFreeJVHost( vcsbeam_context *vm );
void vmFreeJVDevice( vcsbeam_context *vm );

void vmMallocEHost( vcsbeam_context *vm );
void vmMallocEDevice( vcsbeam_context *vm );
void vmFreeEHost( vcsbeam_context *vm );
void vmFreeEDevice( vcsbeam_context *vm );

void vmMallocSHost( vcsbeam_context *vm );
void vmMallocSDevice( vcsbeam_context *vm );
void vmFreeSHost( vcsbeam_context *vm );
void vmFreeSDevice( vcsbeam_context *vm );

void vmMallocJHost( vcsbeam_context *vm );
void vmMallocJDevice( vcsbeam_context *vm );
void vmFreeJHost( vcsbeam_context *vm );
void vmFreeJDevice( vcsbeam_context *vm );

void vmMallocDHost( vcsbeam_context *vm );
void vmMallocDDevice( vcsbeam_context *vm );
void vmFreeDHost( vcsbeam_context *vm );
void vmFreeDDevice( vcsbeam_context *vm );

void vmMallocPQIdxsHost( vcsbeam_context *vm );
void vmMallocPQIdxsDevice( vcsbeam_context *vm );
void vmFreePQIdxsHost( vcsbeam_context *vm );
void vmFreePQIdxsDevice( vcsbeam_context *vm );


void vmSetMaxGPUMem( vcsbeam_context *vm, int nchunks );
void vmPushChunk( vcsbeam_context *vm );
void vmPushJ( vcsbeam_context *vm );

/**
 * vmCreatePrimaryBeam
 * ===================
 *
 * Allocates memory for the primary beam matrices ("B")
 * (see Eq. (30) in Ord et al. (2019))
 */
void vmCreatePrimaryBeam( vcsbeam_context *vm );

/**
 * vmCreateGeometricDelays
 * =======================
 *
 * Allocates memory for the geometric delay arrays ("phi") on both host and device.
 * Free with free_geometric_delays()
 */
void vmCreateGeometricDelays( vcsbeam_context *vm );

void vmCreateCudaStreams( vcsbeam_context *vm );
void vmDestroyCudaStreams( vcsbeam_context *vm );

void vmCreateStatistics( vcsbeam_context *vm, mpi_psrfits *mpfs );
void vmDestroyStatistics( vcsbeam_context *vm );

void vmSetPolIdxLists( vcsbeam_context *vm );
void vmCalcJ( vcsbeam_context *vm );
void vmCalcJonesAndDelays( vcsbeam_context *vm, double *ras_hours, double *decs_degs, beam_geom *beam_geom_vals );

void vmParsePointingFile( vcsbeam_context *vm, const char *filename );
void vmSetNumPointings( vcsbeam_context *vm, unsigned int npointings );

void vmParseFlaggedTilenamesFile( char *filename, calibration *cal ); // (defined in calibration.c)
void vmSetCustomTileFlags( vcsbeam_context *vm ); // (defined in calibration.c);

// (defined in calibration.c:)
void vmLoadRTSSolution( vcsbeam_context *vm );
void vmLoadOffringaSolution( vcsbeam_context *vm );

void vmCreateFilenames( vcsbeam_context *vm );
void vmGetVoltFilename( vcsbeam_context *vm, unsigned int coarse_chan_idx, uint64_t gps_second, char *filename );
void vmGetLegacyVoltFilename( vcsbeam_context *vm, unsigned int coarse_chan_idx, uint64_t gps_second, char *filename );
void vmDestroyFilenames( vcsbeam_context *vm );

void vmLoadObsMetafits( vcsbeam_context *vm, char *filename );
void vmLoadCalMetafits( vcsbeam_context *vm, char *filename );

void vmApplyCalibrationCorrections( vcsbeam_context *vm );

void vmGetVoltageMetadata( vcsbeam_context *vm );

long unsigned int get_relative_gps( MetafitsMetadata *obs_metadata, long int relative_begin );
long unsigned int parse_begin_string( MetafitsMetadata *obs_metadata, char *begin_str );
uintptr_t parse_coarse_chan_string( MetafitsMetadata *obs_metadata, char *begin_coarse_chan_str );

void vmSetNumNotFlaggedRFInputs( vcsbeam_context *vm );

Rfinput *find_matching_rf_input( MetafitsMetadata *metadata, Rfinput *rfinput );
Antenna *find_matching_antenna( MetafitsMetadata *metadata, Rfinput *rfinput );
Antenna *find_antenna_by_name( MetafitsMetadata *obs_metadata, char *tile_name );

void get_mwalib_version( char *version_str );

/**
 * vmLoadFilter
 * ============
 *
 * Load a set of filter coefficients
 * Inputs:
 *   FILTERNAME - string specifying a filter. There should be a corresponding
 *                file in the RUNTIME_DIR called FILTERNAME.dat.
 *   TYPE       - Whether it's an ANALYSIS_FILTER or a SYNTHESIS_FILTER
 *   NCHANS     - the number of channels that this filter will be applied to.
 *                For both ANALYSIS and SYNTHESIS filters, this should be
 *                the number of ANALYSIS channels.
 */
void vmLoadFilter( vcsbeam_context *vm, char *filtername, filter_type type, int nchans );


/**
 * FREE_PFB_FILTER
 * ===============
 *
 * Free the memory allocated in vmLoadFilter()
 */
void free_pfb_filter( pfb_filter *filter );


host_buffer *vmInitReadBuffer( size_t read_size, size_t margin_size );
vm_error vmReadBufferCopyMargin( host_buffer *rb );
void vmFreeReadBuffer( host_buffer *rb );

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

/* In order to communicate easily with the GPU, the "phi" array, which
 * contains the complex-valued geometric delay terms, is implemented
 * as a 1D array, which uses the following indexing macro to access
 * the correct term for a given (p)ointing, (a)ntenna, and (c)hannel.
 */

#define PHI_IDX(p,a,c,na,nc)   ((p) * (nc)*(na)  + \
                                (a) * (nc)       + \
                                (c))

#ifdef __cplusplus
extern "C" {
#endif

/* Calculate the geometric delay (in radians) for the given pointings
 */
void vmCalcPhi(
        vcsbeam_context   *vm,
        beam_geom         *beam_geom_vals );

void calc_geometric_delays(
        beam_geom         *beam_geom_vals,
        uint32_t           freq_hz,
        MetafitsMetadata  *obs_metadata,
        gpuDoubleComplex   *phi );


/* Free memory allocated with create_geometric_delays()
 */
void free_geometric_delays( geometric_delays *gdelays );

/* Copy host memory block to device
 */
void vmPushPhi( vcsbeam_context *vm );

double calc_array_factor(
        MetafitsMetadata *obs_metadata,
        uint32_t          freq_hz,
        beam_geom        *bg1,
        beam_geom        *bg2 );

void calc_beam_geom(
        double            ras_hours,
        double            decs_degs,
        double            mjd,
        beam_geom        *bg );

void dec2hms( char *out, double in, int sflag );
void utc2mjd( char *, double *, double * ); // "2000-01-01T00:00:00" --> MJD_int + MJD_fraction
void mjd2lst( double, double * );

double parse_dec( char* ); // "01:23:45.67" --> Dec in degrees
double parse_ra( char* );  // "01:23:45.67" --> RA  in degrees

void vmScaleFilterCoeffs( vcsbeam_context *vm, filter_type type, double scale_factor );

void vmInitForwardPFB( vcsbeam_context *vm, int M, pfb_flags flags );
void vmUploadForwardPFBChunk( vcsbeam_context *vm );
void vmWOLAChunk( vcsbeam_context *vm );
void vmFPGARoundingChunk( vcsbeam_context *vm );
void vmFFTChunk( vcsbeam_context *vm );
void vmPackChunk( vcsbeam_context *vm );
void vmDownloadForwardPFBChunk( vcsbeam_context *vm );
void vmExecuteForwardPFB( vcsbeam_context *vm );

void vmWritePFBOutputToFile( vcsbeam_context *vm );

void vmFreeForwardPFB( forward_pfb *fpfb );

void vmReportPerformanceStats( vcsbeam_context *vm );

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

void cu_invert_pfb( gpuDoubleComplex *data_buffer_fine, int file_no,
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

#define REAL_NIBBLE_TO_UINT8(X)  (((X) >> 4) & 0xf)
#define IMAG_NIBBLE_TO_UINT8(X)  ((X) & 0xf)
#define UINT8_TO_INT(X)          ((X) >= 0x8 ? (signed int)(X) - 0x10 : (signed int)(X))
#define RE_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))))
#define IM_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X))))
#define UCMPLX4_TO_CMPLX_FLT(X)  (make_gpuDoubleComplex((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))), \
                                         (float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X)))))
#define DETECT(X)                (gpuCreal(gpuCmul(X,gpuConj(X))))


#ifdef __cplusplus
extern "C" {
#endif


void cu_form_incoh_beam(
        uint8_t *data, uint8_t *d_data, size_t data_size,
        float *d_incoh,
        unsigned int nsample, int nchan, int ninput,
        float *offsets, float *d_offsets,
        float *scales, float *d_scales,
        uint8_t *Iscaled, uint8_t *d_Iscaled, size_t Iscaled_size
        );


void vmApplyJChunk( vcsbeam_context *vm );
void vmBeamformChunk( vcsbeam_context *vm );
void vmBeamformSecond( vcsbeam_context *vm );
void vmPullE( vcsbeam_context *vm );
void vmPullS( vcsbeam_context *vm );

void prepare_data_buffer_fine( gpuDoubleComplex *data_buffer_fine, vcsbeam_context *vm,
                    uintptr_t timestep_idx );
void vmSendSToFits( vcsbeam_context *vm, mpi_psrfits *mpfs );

float *create_pinned_data_buffer( size_t size );

void vmPushPolIdxLists( vcsbeam_context *vm );

vm_error vmReadNextSecond( vcsbeam_context *vm );

gpuDoubleComplex ****create_invJi( int nstation, int nchan, int pol );
void              destroy_invJi( gpuDoubleComplex ****array, int nstation, int nchan, int npol );

gpuDoubleComplex *create_data_buffer_fine( int npointing, int nsamples, int nchan, int npol );

void allocate_input_output_arrays( void **data, void **d_data, size_t size );
void free_input_output_arrays( void *data, void *d_data );

#ifdef __cplusplus
}
#endif



/********************
 *                  *
 *   BEAM_PSRFITS   *
 *                  *
 ********************/

#ifdef __cplusplus
extern "C" {
#endif

void populate_psrfits_header(
        vcsbeam_context  *vm,
        struct psrfits   *pf,
        int               max_sec_per_file,
        int               outpol,
        beam_geom        *beam_geom_vals,
        char             *incoh_basename,
        bool              is_coherent );

void populate_spliced_psrfits_header(
        vcsbeam_context  *vm,
        struct psrfits   *pf,
        int               max_sec_per_file,
        int               outpol,
        beam_geom        *beam_geom_vals,
        char             *basename,
        bool              is_coherent );

void free_psrfits( struct psrfits *pf );

void vmInitMPIPsrfits(
        vcsbeam_context *vm,
        mpi_psrfits *mpf,
        int max_sec_per_file,
        int nstokes,
        beam_geom *bg,
        char *outfile,
        bool is_coherent );

void free_mpi_psrfits( mpi_psrfits *mpf );

void gather_splice_psrfits( mpi_psrfits *mpf );
void wait_splice_psrfits( mpi_psrfits *mpf );

#ifdef __cplusplus
}
#endif


/*******************
 *                 *
 *    BEAM_VDIF    *
 *                 *
 *******************/

/* convenience type - this just collects all the vdif info together */
struct vdifinfo
{
    int frame_length;         // Length of the vdif frame
    int frame_rate;           // Frames per second
    size_t samples_per_frame; // Number of time samples per vdif frame
    int dataarraylength;      // The frame length minus the header
    int bits;                 // Bits per sample
    int nchan;                // Channels per frame
    int chan_width;           // Channel width in hertz
    int sample_rate;          // Sample rate in hertz
    int npol;                 // Number of polarisations
    int iscomplex;            // Complex sample flag
    int threadid;             // Which thread are we
    char stationid[3];        // Which station are we
    char exp_name[17];        // Experiment name
    char scan_name[17];       // Scan_name
    size_t block_size;        // Size of one second of output data including headers
    size_t sizeof_buffer;     // Size of one second of 32bit complex beam data (no headers)
    size_t sizeof_beam;       // Size of 1 sample of 32bit complex beam data (no headers)
    float *b_scales;          // Bandpass mean
    float *b_offsets;         // Bandpass offset

    // observation info
    char telescope[24];
    char source[24];
    char obs_mode[8];
    double fctr;

    char ra_str[16];
    char dec_str[16];

    double BW;

    double MJD_start;
    double sec_offset;
    double MJD_epoch;

    char date_obs[24];
    char basefilename[1024];

    int got_scales;
};

#ifdef __cplusplus
extern "C" {
#endif


void vdif_write_data( struct vdifinfo *vf, int8_t *output );
void vdif_write_second( struct vdifinfo *vf, vdif_header *vhdr,
        float *data_buffer_vdif );

void vmPopulateVDIFHeader(
        vcsbeam_context  *vm,
        beam_geom        *beam_geom_vals,
        double           mjd_start,
        double           sec_offset );

void to_offset_binary( int8_t *i, int n );


#ifdef __cplusplus
}
#endif



/********************
 *                  *
 *   PRIMARY_BEAM   *
 *                  *
 ********************/

#define NCONFIGS           138
#define DEAD_CONFIG        (NCONFIGS - 1)
#define MANY_DEAD_DIPOLES  -1

#define PB_IDX(p,a,pol,na,npol) ((p)*(na)*(npol) + (a)*(npol) + (pol))


#ifdef __cplusplus
extern "C" {
#endif


void handle_hyperbeam_error(char file[], int line_num, const char function_name[]);

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

void vmCalcB(
        vcsbeam_context   *vm,
        beam_geom         *beam_geom_vals );

void free_primary_beam( primary_beam *pb );

void calc_normalised_beam_response( FEEBeam *beam, double az, double za, double freq_hz, uint32_t *delays, double *amps, double *IQUV, gpuDoubleComplex *J, bool apply_pa_correction );

#ifdef __cplusplus
}
#endif




#endif
