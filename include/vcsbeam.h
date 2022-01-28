/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __VCSBEAM_H__
#define __VCSBEAM_H__


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cufft.h>

#include <mwalib.h>

#cmakedefine VCSBEAM_VERSION  "@VCSBEAM_VERSION@"
#cmakedefine RUNTIME_DIR      "@RUNTIME_DIR@"

#define ERROR_MESSAGE_LEN  1024
#define MAX_COMMAND_LENGTH 1024

// Calibration solution types
#define CAL_NONE      0
#define CAL_RTS       1
#define CAL_OFFRINGA  2



#define PERFORMANCE_MAX_NUM_STOPWATCHES  16
#define PERFORMANCE_MAX_START_STOP       4096



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
    cuDoubleComplex *twiddles; // twiddle factors
    filter_type      type;
} pfb_filter;

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
    cuDoubleComplex   *phi;
    cuDoubleComplex   *d_phi;
    uintptr_t          npointings;
    uintptr_t          nant;
    uintptr_t          nchan;
    double            *chan_freqs_hz;
    MetafitsMetadata  *obs_metadata;
} geometric_delays;

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

    cuFloatComplex   *d_weighted_overlap_add; // A "temporary" array on the device for mid-calculation product
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
    cufftHandle       plan;                   // The cuFFT plan for performing the FFT part of the forward PFB
} forward_pfb;


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



/* TODO: Move this (and its pointer in vcsbeam_context) out of the public header */
struct vdifinfo
{
    int frame_length;         // Length of the vdif frame
    int frame_rate;           // Frames per second
    size_t samples_per_frame; // Number of time samples per vdif frame
    int dataarraylength;      // The frame length minus the header
    int bits;                 // Bits per sample
    int nchan;                // Channels per framme
    int chan_width;           // Channel width in hertz
    int sample_rate;          // Sample rate in hertz
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

    long double MJD_epoch;

    char date_obs[24];
    char basefilename[1024];

    int got_scales;
};



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

    cuDoubleComplex *e, *d_e;         // The buffers for the beamformed voltages on host/device
    uintptr_t e_size_bytes;           // The size of S in bytes
    uintptr_t d_e_size_bytes;         // The size of d_S in bytes

    cuDoubleComplex *J, *d_J;         // The buffers for Jones matrices on host/device
    uintptr_t J_size_bytes;           // The size of J in bytes
    uintptr_t d_J_size_bytes;         // The size of d_J in bytes

    cuDoubleComplex *Jv_P, *d_Jv_P;   // The buffers for Jones-corrected voltages on host/device
    cuDoubleComplex *Jv_Q, *d_Jv_Q;   // The buffers for Jones-corrected voltages on host/device
    uintptr_t Jv_size_bytes;          // The size of Jv in bytes
    uintptr_t d_Jv_size_bytes;        // The size of d_Jv in bytes

    cuDoubleComplex *D, *d_D;         // The buffers for calibration solutions on host/device
    uintptr_t D_size_bytes;           // The size of D in bytes
    uintptr_t d_D_size_bytes;         // The size of d_D in bytes

    uint32_t *polP_idxs, *d_polP_idxs; // List of indices for VCS-ordered data
    uint32_t *polQ_idxs, *d_polQ_idxs;
    uintptr_t pol_idxs_size_bytes;     // The size of (each of) the P/Q idxs (host) arrays
    uintptr_t d_pol_idxs_size_bytes;   // The size of (each of) the P/Q idxs (device) arrays

    geometric_delays gdelays;         // Geometric delays due to tile layout

    unsigned int npointing;           // Number of requested tied array beam pointings

    uintptr_t max_gpu_mem_bytes;      // The maximum allowed GPU memory to use (in bytes)
    uint32_t chunks_per_second;       // The number of chunks to process on device per second of data
                                      // (data_size_bytes = d_data_size_bytes * chunks_per_second)

    int num_coarse_chans_to_process;  // The number of coarse channels to be processed
    int *coarse_chan_idxs_to_process; // A list of the coarse chan idxs to be processed

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

    // VDIF output
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

    cudaStream_t *streams;            // CUDA streams used for forming multiple tied-array beams

    logger *log;                      // Used for log messages
    char log_message[MAX_COMMAND_LENGTH];

    char error_message[ERROR_MESSAGE_LEN]; // MWALIB error message buffer
} vcsbeam_context;

#ifdef __cplusplus
extern "C" {
#endif

vcsbeam_context *vmInit( bool use_mpi );

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
void destroy_vcsbeam_context( vcsbeam_context *vm );
void vmSetOutputChannelisation( vcsbeam_context *vm, bool out_fine, bool out_coarse );

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

void vmSetMaxGPUMem( vcsbeam_context *vm, uintptr_t max_gpu_mem_bytes );
void vmPushChunk( vcsbeam_context *vm );
void vmPushJ( vcsbeam_context *vm );

void vmCreateGeometricDelays( vcsbeam_context *vm );

void vmCreateCudaStreams( vcsbeam_context *vm );
void vmDestroyCudaStreams( vcsbeam_context *vm );

void vmDestroyStatistics( vcsbeam_context *vm );

void vmSetPolIdxLists( vcsbeam_context *vm );

void vmParsePointingFile( vcsbeam_context *vm, const char *filename );
void vmSetNumPointings( vcsbeam_context *vm, unsigned int npointings );

void vmParseFlaggedTilenamesFile( char *filename, calibration *cal ); // (defined in calibration.c)
void vmSetCustomTileFlags( vcsbeam_context *vm ); // (defined in calibration.c);

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

void vmLoadFilter( vcsbeam_context *vm, char *filtername, filter_type type, int nchans );

host_buffer *vmInitReadBuffer( size_t read_size, size_t margin_size );
vm_error vmReadBufferCopyMargin( host_buffer *rb );
void vmFreeReadBuffer( host_buffer *rb );

void vmCalcPhi(
        vcsbeam_context   *vm,
        beam_geom         *beam_geom_vals );

/* Copy host memory block to device
 */
void vmPushPhi( vcsbeam_context *vm );

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

void vmApplyJChunk( vcsbeam_context *vm );
void vmBeamformChunk( vcsbeam_context *vm );
void vmBeamformSecond( vcsbeam_context *vm );
void vmPullE( vcsbeam_context *vm );
void vmPullS( vcsbeam_context *vm );

void vmPushPolIdxLists( vcsbeam_context *vm );

vm_error vmReadNextSecond( vcsbeam_context *vm );


#ifdef __cplusplus
}
#endif




#endif
