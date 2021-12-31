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

#include "vcsbeam.h"

/**
 * Initialises a VCSBeam context struct.
 *
 * @param use_mpi Set up the struct for using MPI
 * @return A pointer to a newly allocated `vcsbeam_context` struct
 *
 * Once the VCSBeam context is finished with, it should be freed with a call
 * to destroy_vcsbeam_context().
 */
vcsbeam_context *vmInit( bool use_mpi )
{
    // Allocate memory for the VCSBEAM_METADATA struct
    vcsbeam_context *vm = (vcsbeam_context *)malloc( sizeof(vcsbeam_context) );

    // Initialise MPI
    vm->use_mpi = use_mpi;
    if (use_mpi)
    {
        MPI_Init( NULL, NULL );
        MPI_Comm_size( MPI_COMM_WORLD, &vm->mpi_size );
        MPI_Comm_rank( MPI_COMM_WORLD, &vm->mpi_rank );
    }
    else
    {
        vm->mpi_size = 1;
        vm->mpi_rank = PERFORMANCE_NO_MPI;
    }
    vm->writer = 0;

    // TODO: Change this to give user flexibility of how to use mpi structure
    vm->ncoarse_chans   = vm->mpi_size;
    vm->coarse_chan_idx = (vm->use_mpi ? vm->mpi_rank : 0);

    // Start with the first chunk
    vm->chunk_to_load = 0;

    // Initialise data pointers to NULL
    vm->v           = NULL;
    vm->d_v         = NULL;
    vm->S           = NULL;
    vm->d_S         = NULL;
    vm->e           = NULL;
    vm->d_e         = NULL;
    vm->J           = NULL;
    vm->d_J         = NULL;
    vm->Jv_P        = NULL;
    vm->d_Jv_P      = NULL;
    vm->Jv_Q        = NULL;
    vm->d_Jv_Q      = NULL;
    vm->D           = NULL;
    vm->d_D         = NULL;
    vm->polP_idxs   = NULL;
    vm->polQ_idxs   = NULL;
    vm->d_polP_idxs = NULL;
    vm->d_polQ_idxs = NULL;

    vm->d_v_size_bytes        = 0;
    vm->pol_idxs_size_bytes   = 0;
    vm->d_pol_idxs_size_bytes = 0;
    vm->D_size_bytes          = 0;
    vm->d_D_size_bytes        = 0;
    vm->J_size_bytes          = 0;
    vm->d_J_size_bytes        = 0;
    vm->e_size_bytes          = 0;
    vm->d_e_size_bytes        = 0;
    vm->S_size_bytes          = 0;
    vm->d_S_size_bytes        = 0;
    vm->Jv_size_bytes         = 0;
    vm->d_Jv_size_bytes       = 0;

    // Default: data is legacy VCS format (VM_INT4)
    vm->datatype = VM_INT4;

    // Calibration
    init_calibration( &vm->cal );

    // Initially not bound to any observations
    vm->obs_context  = NULL;
    vm->obs_metadata = NULL;

    vm->vcs_context  = NULL;
    vm->vcs_metadata = NULL;

    vm->cal_context  = NULL;
    vm->cal_metadata = NULL;

    // Initially, do nothing
    vm->do_forward_pfb = false;
    vm->do_inverse_pfb = false;

    vm->output_fine_channels = false;
    vm->output_coarse_channels = false;

    // No filters
    vm->analysis_filter = NULL;
    vm->synth_filter    = NULL;

    // No forward PFB
    vm->fpfb = NULL;

    // Start (reading) at the beginning
    vm->current_gps_idx = 0;

    // No input
    vm->datadir   = NULL;
    vm->filenames = NULL;
    vm->nfiles    = 0;

    // Start a logger
    vm->log = create_logger( stdout, vm->mpi_rank );
    logger_add_stopwatch( vm->log, "read",      "Reading in data" );
    logger_add_stopwatch( vm->log, "upload",    "Uploading the data to the device" );
    logger_add_stopwatch( vm->log, "pfb",       "Performing the PFB" );
    logger_add_stopwatch( vm->log, "pfb-wola",  "Weighted overlap-add" );
    logger_add_stopwatch( vm->log, "pfb-round", "FPGA rounding and demotion" );
    logger_add_stopwatch( vm->log, "pfb-fft",   "Performing FFT" );
    logger_add_stopwatch( vm->log, "pfb-pack",  "Packing the data into the recombined format" );
    logger_add_stopwatch( vm->log, "delay",     "Calculating geometric and cable delays" );
    logger_add_stopwatch( vm->log, "calc",      "Calculating tied-array beam" );
    logger_add_stopwatch( vm->log, "ipfb",      "Inverting the PFB" );
    logger_add_stopwatch( vm->log, "download",  "Downloading the data to the host" );
    logger_add_stopwatch( vm->log, "splice",    "Splicing coarse channels together" );
    logger_add_stopwatch( vm->log, "write",     "Writing out data to file" );

    // Initialise pointing RAs and Decs to NULL
    vm->ras_hours = NULL;
    vm->decs_degs = NULL;

    // Return the new struct pointer
    return vm;
}

/**
 * Binds a set of observation files to the VCSBeam context.
 *
 * @param first_coarse_chan_str A string representation of the lowest coarse
 *        channel to be processed
 * @param num_coarse_chans_to_process The number of coarse channels to be
 *        processed
 * @param coarse_chan_idx_offset The first coarse channel (relative to
 *        `first_coarse_chan_str`) to be processed by this MPI process
 * @param first_gps_second_str A string representation of the first GPS second
 *        to be processed
 * @param num_gps_seconds_to_process The number of GPS seconds to be processed
 * @param gps_second_offset The first GPS second (relative to
 *        `first_gps_second_str`) to be processed by this MPI process
 * @param datadir The directory containing the observation data files
 */
void vmBindObsData(
        vcsbeam_context *vm,
        char *first_coarse_chan_str, int num_coarse_chans_to_process, int coarse_chan_idx_offset,
        char *first_gps_second_str, int num_gps_seconds_to_process, int gps_second_offset,
        char *datadir )
{
    // Convert the (first) chan and gps strings to numbers
    uint32_t first_gps_second = parse_begin_string( vm->obs_metadata, first_gps_second_str ) + gps_second_offset;
    int first_coarse_chan_idx = parse_coarse_chan_string( vm->obs_metadata, first_coarse_chan_str ) + coarse_chan_idx_offset;

    // Get the voltage context and metadata
    vm->num_gps_seconds_to_process  = num_gps_seconds_to_process;
    vm->num_coarse_chans_to_process = num_coarse_chans_to_process;

    // Construct the gps second and coarse chan idx arrays
    vm->gps_seconds_to_process = (uint32_t *)malloc( vm->num_gps_seconds_to_process * sizeof(uint32_t) );
    vm->coarse_chan_idxs_to_process = (int *)malloc( vm->num_coarse_chans_to_process * sizeof(int) );

    int g;
    for (g = 0; g < vm->num_gps_seconds_to_process; g++)
    {
        vm->gps_seconds_to_process[g] = g + first_gps_second;
    }

    int c;
    for (c = 0; c < num_coarse_chans_to_process; c++)
    {
        vm->coarse_chan_idxs_to_process[c] = c + first_coarse_chan_idx;
    }

    // Copy across the data directory
    vm->datadir = (char *)malloc( strlen( datadir ) + 1 );
    strcpy( vm->datadir, datadir );

    vmGetVoltageMetadata( vm );

    // Set the default output to match the same channelisation as the input
    vmSetOutputChannelisation( vm,
            (vm->obs_metadata->mwa_version == VCSLegacyRecombined),
            (vm->obs_metadata->mwa_version == VCSMWAXv2)
            );

    // Compute the "shorthand" variables
    vm->sample_rate = vm->vcs_metadata->num_samples_per_voltage_block *
                      vm->vcs_metadata->num_voltage_blocks_per_second;

    vm->bytes_per_second = vm->vcs_metadata->num_voltage_blocks_per_second *
                           vm->vcs_metadata->voltage_block_size_bytes;

    vm->nchan = vm->obs_metadata->num_volt_fine_chans_per_coarse;

    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
    {
        vm->fine_sample_rate = vm->sample_rate;
        vm->nfine_chan       = vm->nchan;
    }
    else
    {
        // (These will get set if doing a forward PFB)
        vm->fine_sample_rate = 0;
        vm->nfine_chan       = 0;
    }

    // Actually allocate the read buffer
    vmMallocVHost( vm );
}

/**
 * Binds a calibration solution to the VCSBeam context.
 *
 * @param caldir The directory containing RTS solution files, OR the path of
 *        an Offringa-style calibration solution file
 * @param cal_type Either `CAL_RTS` or `CAL_OFFRINGA`
 * @param use_bandpass Whether to include the Bandpass information (relevant
 *        for RTS solutions only)
 * @param flags_file A file containing names of (extra) tiles to be flagged
 */
void vmBindCalibrationData( vcsbeam_context *vm,
        char   *caldir,
        int     cal_type,
        bool    use_bandpass,  // only relevant for RTS observations
        char   *flags_file )
{
    // Set the calibration type (either RTS or Offringa)
    vm->cal.cal_type     = cal_type;
    vm->cal.use_bandpass = use_bandpass;

    vm->cal.caldir       = (char *)malloc( strlen( caldir ) + 1 );
    strcpy( vm->cal.caldir, caldir );

    vm->cal.flags_file   = (char *)malloc( strlen( flags_file ) + 1 );
    strcpy( vm->cal.flags_file, flags_file );
}

/**
 * Reads in a calibration solution.
 *
 * Calls either vmLoadRTSSolution() or vmLoadOffringaSolution() depending
 * on whether `vm&rarr;cal.cal_type` is set to `CAL_RTS` or `CAL_OFFRINGA`.
 * Afterwards, vmSetCustomTileFlags() is called.
 */
void vmReadCalibration( vcsbeam_context *vm )
{
    // Read in the calibration data from file
    if (vm->cal.cal_type == CAL_RTS)
        vmLoadRTSSolution( vm );
    else if (vm->cal.cal_type == CAL_OFFRINGA)
        vmLoadOffringaSolution( vm );

    // Flag extra tiles that need flagging
    vmSetCustomTileFlags( vm );
}

/**
 * Finds an antenna in an observation with a given name.
 *
 * @param obs_metadata The observation metadata to be searched
 * @param tile_name The name of the tile being sought
 * @return A pointer to an mwalib Antenna struct with a matching tile name
 *
 * If no tile with the given name is found in the observation, `NULL` is
 * returned.
 */
Antenna *find_antenna_by_name( MetafitsMetadata *obs_metadata, char *tile_name )
{
    int i;
    for (i = 0; i < (int)(obs_metadata->num_ants); i++)
    {
        if (strcmp( tile_name, obs_metadata->antennas[i].tile_name ) == 0)
            return &(obs_metadata->antennas[i]);
    }

    // If we got this far, no antenna with that name was found
    return NULL;
}

/**
 * Frees the memory associated with the VCSBeam context.
 *
 * After freeing the memory associated with the VCSBeam context's member
 * variables, this function frees the VCSBeam context itself.
 *
 * @todo Rename this function to a `vm...` name
 */
void destroy_vcsbeam_context( vcsbeam_context *vm )
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

    destroy_logger( vm->log );

    // Free the RA and Dec pointing arrays
    if (vm->ras_hours != NULL)  free( vm->ras_hours );
    if (vm->decs_degs != NULL)  free( vm->decs_degs );

    // Calibration
    free_calibration( &vm->cal );

    // Filters
    if (vm->analysis_filter != NULL)
        free_pfb_filter( vm->analysis_filter );

    if (vm->synth_filter != NULL)
        free_pfb_filter( vm->synth_filter );

    // Forward PFB
    if (vm->fpfb != NULL)
        vmFreeForwardPFB( vm->fpfb );

    // Read buffer
    vmFreeVHost( vm );

    // Input data info
    if (vm->datadir != NULL)
        free( vm->datadir );

    vmDestroyFilenames( vm );

    // Finalise MPI
    if (vm->use_mpi)
        MPI_Finalize();

    // Finally, free the struct itself
    free( vm );
}

/**
 * Sets flags governing whether the PFB and inverse PFB routines are run
 * depending on the input and output channelisations.
 *
 * @param out_fine Sets the flag for fine channelised output
 * @param out_coarse Sets the flag for coarse_channelised output
 *
 * Whether the (forward) PFB or the inverse PFB needs to be run depends on the
 * input channelisation (fine = Legacy, or coarse = MWAX), and what
 * channelisations are desired in output (fine or coarse).
 * The `vm&rarr;do_forward_pfb` and `vm&rarr;do_inverse_pfb` are set
 * accordingly.
 *
 * The following table describes all possible scenarios:
 * | Input mode | Output fine? | Output coarse? | Do forward PFB? | Do inverse PFB? |
 * | :--------: | :----------: | :------------: | :-------------: | :-------------: |
 * | Legacy     | no           | no             | no              | no              |
 * | Legacy     | no           | yes            | no              | yes             |
 * | Legacy     | yes          | no             | no              | no              |
 * | Legacy     | yes          | yes            | no              | yes             |
 * | MWAX       | no           | no             | no              | no              |
 * | MWAX       | no           | yes            | yes             | yes             |
 * | MWAX       | yes          | no             | yes             | no              |
 * | MWAX       | yes          | yes            | yes             | yes             |
 *
 * Note that both forward and inverse PFBs are required for MWAX data even
 * when fine-channelised output is not requested.
 * This is because the beamforming operation requires sufficiently fine
 * channels in order to avoid decoherence across the channels.
 */
void vmSetOutputChannelisation( vcsbeam_context *vm, bool out_fine, bool out_coarse )
{
    vm->output_fine_channels   = out_fine;
    vm->output_coarse_channels = out_coarse;

    if (!out_fine && !out_coarse) // #1 and #5
    {
        vm->do_forward_pfb = false;
        vm->do_inverse_pfb = false;
        return;
    }

    // The inverse pfb is needed iff coarse channels are output
    vm->do_inverse_pfb = vm->output_coarse_channels;

    // The forward pfb is only needed for the (remaining) MWAX options
    vm->do_forward_pfb = (vm->obs_metadata->mwa_version == VCSMWAXv2);
}

void vmMallocVHost( vcsbeam_context *vm )
{
    // If the input is legacy, the read buffer can be just enough
    // for one second's worth of data.
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
        vm->v = vmInitReadBuffer( vm->bytes_per_second, 0 );
    // However, for MWAX, the read buffer must be slightly bigger
    // to accommodate the PFB's extra taps.
    // For now, this is FIXED to one voltage block. Changing this
    // will break the gpu kernels that do the fine PFB, for which
    // the offset into the data is currently hard-coded.
    else if (vm->obs_metadata->mwa_version == VCSMWAXv2)
        vm->v = vmInitReadBuffer( vm->bytes_per_second, vm->vcs_metadata->voltage_block_size_bytes );
}

void vmFreeVHost( vcsbeam_context *vm )
{
    vmFreeReadBuffer( vm->v );
}

/**
 * vmMallocVDevice
 * ===============
 *
 * Allocates device memory for the input voltages for the beamformer.
 * If the observation is legacy, then allocates new memory.
 * If the observation is MWAX, then just point to the assumed
 * already existing memory allocated in VM->FPFB.
 */
void vmMallocVDevice( vcsbeam_context *vm )
{
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
    {
        vm->d_v_size_bytes = vm->bytes_per_second / vm->chunks_per_second;
        cudaMalloc( (void **)&vm->d_v,  vm->d_v_size_bytes );
        cudaCheckErrors( "vmMallocVDevice: cudaMalloc failed" );
    }
    else // if (vm->obs_metadata->mwa_version == VCSMWAXv2)
    {
        vm->d_v            = vm->fpfb->d_vcs_data;
        vm->d_v_size_bytes = vm->fpfb->d_vcs_size;
        // TODO: ^^^ Check that these are already the same thing
    }
}

/**
 * vmFreeVDevice
 * =============
 *
 * Frees the memory allocated with vmMallocVDevice.
 * If the observation is legacy, then free the memory.
 * If the observation is MWAX, then do nothing; the memory
 * should be freed via a call to vmFreeForwardPFB.
 */
void vmFreeVDevice( vcsbeam_context *vm )
{
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
    {
        cudaFree( vm->d_v );
        cudaCheckErrors( "vmFreeVDevice: cudaFree failed" );
    }
}

void vmMallocJVHost( vcsbeam_context *vm )
{
    vm->Jv_size_bytes =
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->fine_sample_rate *
        sizeof(cuDoubleComplex);

    cudaMallocHost( (void **)&(vm->Jv_P), vm->Jv_size_bytes );
    cudaCheckErrors( "vmMallocJVHost: cudaMallocHost(Jv_P) failed" );
    cudaMallocHost( (void **)&(vm->Jv_Q), vm->Jv_size_bytes );
    cudaCheckErrors( "vmMallocJVHost: cudaMallocHost(Jv_Q) failed" );
}

void vmFreeJVHost( vcsbeam_context *vm )
{
    cudaFreeHost( vm->Jv_P );
    cudaCheckErrors( "vmFreeJVHost: cudaFreeHost(Jv_P) failed" );
    cudaFreeHost( vm->Jv_Q );
    cudaCheckErrors( "vmFreeJVHost: cudaFreeHost(Jv_Q) failed" );
}

void vmMallocJVDevice( vcsbeam_context *vm )
{
    vm->d_Jv_size_bytes = 
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->fine_sample_rate *
        sizeof(cuDoubleComplex) /
        vm->chunks_per_second;

    cudaMalloc( (void **)&vm->d_Jv_P,  vm->d_Jv_size_bytes );
    cudaCheckErrors( "vmMallocJVDevice: cudaMalloc(d_Jv_P) failed" );
    cudaMalloc( (void **)&vm->d_Jv_Q,  vm->d_Jv_size_bytes );
    cudaCheckErrors( "vmMallocJVDevice: cudaMalloc(d_Jv_Q) failed" );
}

void vmFreeJVDevice( vcsbeam_context *vm )
{
    cudaFree( vm->d_Jv_P );
    cudaCheckErrors( "vmFreeJVDevice: cudaFree(d_Jv_P) failed" );
    cudaFree( vm->d_Jv_Q );
    cudaCheckErrors( "vmFreeJVDevice: cudaFree(d_Jv_Q) failed" );
}

void vmMallocEHost( vcsbeam_context *vm )
{
    vm->e_size_bytes =
        vm->npointing *
        vm->fine_sample_rate *
        vm->nfine_chan *
        vm->obs_metadata->num_ant_pols *
        sizeof(cuDoubleComplex);

    // Allocate memory on host
    cudaMallocHost( (void **)&(vm->e), vm->e_size_bytes );
    cudaCheckErrors( "vmMallocEHost: cudaMallocHost failed" );
}

void vmFreeEHost( vcsbeam_context *vm )
{
    cudaFreeHost( vm->e );
    cudaCheckErrors( "vmFreeEHost: cudaFreeHost failed" );
}

void vmMallocEDevice( vcsbeam_context *vm )
{
    vm->d_e_size_bytes =
        vm->npointing *
        vm->fine_sample_rate *
        vm->nfine_chan *
        vm->obs_metadata->num_ant_pols *
        sizeof(cuDoubleComplex);

    // Allocate memory on device
    cudaMalloc( (void **)&(vm->d_e), vm->d_e_size_bytes );
    cudaCheckErrors( "vmMallocEDevice: cudaMalloc failed" );
}

void vmFreeEDevice( vcsbeam_context *vm )
{
    cudaFree( vm->d_e );
    cudaCheckErrors( "vmFreeEDevice: cudaFree failed" );
}

void vmMallocSHost( vcsbeam_context *vm )
{
    vm->S_size_bytes = vm->npointing * vm->nfine_chan * NSTOKES * vm->fine_sample_rate * sizeof(float);

    // Allocate memory on host
    cudaMallocHost( (void **)&(vm->S), vm->S_size_bytes );
    cudaCheckErrors( "vmMallocSHost: cudaMallocHost failed" );
}

void vmFreeSHost( vcsbeam_context *vm )
{
    cudaFreeHost( vm->S );
    cudaCheckErrors( "vmFreeSHost: cudaFreeHost failed" );
}

void vmMallocSDevice( vcsbeam_context *vm )
{
    vm->d_S_size_bytes = vm->npointing * vm->nfine_chan * NSTOKES * vm->fine_sample_rate * sizeof(float);

    // Allocate memory on device
    cudaMalloc( (void **)&(vm->d_S), vm->d_S_size_bytes );
    cudaCheckErrors( "vmMallocSDevice: cudaMalloc failed" );
}

void vmFreeSDevice( vcsbeam_context *vm )
{
    cudaFree( vm->d_S );
    cudaCheckErrors( "vmFreeSDevice: cudaFree failed" );
}

void vmMallocJHost( vcsbeam_context *vm )
{

    vm->J_size_bytes =
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(cuDoubleComplex);

    // Allocate memory on device
    cudaMallocHost( (void **)&(vm->J), vm->J_size_bytes );
    cudaCheckErrors( "vmMallocJHost: cudaMallocHost failed" );
}

void vmFreeJHost( vcsbeam_context *vm )
{
    cudaFreeHost( vm->J );
    cudaCheckErrors( "vmFreeJHost: cudaFreeHost failed" );
}

void vmMallocJDevice( vcsbeam_context *vm )
{
    vm->d_J_size_bytes =
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(cuDoubleComplex);

    // Allocate memory on device
    cudaMalloc( (void **)&(vm->d_J), vm->d_J_size_bytes );
    cudaCheckErrors( "vmMallocJDevice: cudaMalloc failed" );
}

void vmFreeJDevice( vcsbeam_context *vm )
{
    cudaFree( vm->d_J );
    cudaCheckErrors( "vmFreeJDevice: cudaFree failed" );
}


void vmMallocDHost( vcsbeam_context *vm )
{
    vm->D_size_bytes =
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(cuDoubleComplex);

    // Allocate memory on device
    cudaMallocHost( (void **)&(vm->D), vm->D_size_bytes );
    cudaCheckErrors( "vmMallocDHost: cudaMallocHost failed" );
}

void vmFreeDHost( vcsbeam_context *vm )
{
    cudaFreeHost( vm->D );
    cudaCheckErrors( "vmFreeDHost: cudaFreeHost failed" );
}

void vmMallocDDevice( vcsbeam_context *vm )
{
    vm->d_D_size_bytes =
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(cuDoubleComplex);

    // Allocate memory on device
    cudaMalloc( (void **)&(vm->d_D), vm->d_D_size_bytes );
    cudaCheckErrors( "vmMallocDDevice: cudaMalloc failed" );
}

void vmFreeDDevice( vcsbeam_context *vm )
{
    cudaFree( vm->d_D );
    cudaCheckErrors( "vmFreeDDevice: cudaFree failed" );
}


void vmMallocPQIdxsHost( vcsbeam_context *vm )
{
    vm->pol_idxs_size_bytes = vm->obs_metadata->num_ants * sizeof(uint32_t);

    // Allocate memory on device
    cudaMallocHost( (void **)&(vm->polP_idxs), vm->pol_idxs_size_bytes );
    cudaCheckErrors( "vmMallocPQIdxsHost: cudaMallocHost(polP_idxs) failed" );
    cudaMallocHost( (void **)&(vm->polQ_idxs), vm->pol_idxs_size_bytes );
    cudaCheckErrors( "vmMallocPQIdxsHost: cudaMallocHost(polQ_idxs) failed" );
}

void vmFreePQIdxsHost( vcsbeam_context *vm )
{
    cudaFreeHost( vm->polP_idxs );
    cudaCheckErrors( "vmFreePQIdxsHost: cudaFreeHost(polP_idxs) failed" );
    cudaFreeHost( vm->polQ_idxs );
    cudaCheckErrors( "vmFreePQIdxsHost: cudaFreeHost(polQ_idxs) failed" );
}


void vmMallocPQIdxsDevice( vcsbeam_context *vm )
{
    vm->d_pol_idxs_size_bytes = vm->obs_metadata->num_ants * sizeof(uint32_t);

    // Allocate memory on device
    cudaMalloc( (void **)&(vm->d_polP_idxs), vm->d_pol_idxs_size_bytes );
    cudaCheckErrors( "vmMallocPQIdxsDevice: cudaMalloc(polP_idxs) failed" );
    cudaMalloc( (void **)&(vm->d_polQ_idxs), vm->d_pol_idxs_size_bytes );
    cudaCheckErrors( "vmMallocPQIdxsDevice: cudaMalloc(polQ_idxs) failed" );
}

void vmFreePQIdxsDevice( vcsbeam_context *vm )
{
    cudaFree( vm->d_polP_idxs );
    cudaCheckErrors( "vmFreePQIdxsDevice: cudaFree(polP_idxs) failed" );
    cudaFree( vm->d_polQ_idxs );
    cudaCheckErrors( "vmFreePQIdxsDevice: cudaFree(polQ_idxs) failed" );
}

/**
 * vmSetMaxGPUMem
 * ==============
 *
 * Tries to cleverly figure out how many chunks are needed to
 * fit everything on the GPU.
 *
 * LOGIC IS CURRENTLY FAULTY AND INCOMPLETE. DO NOT USE
 */
void vmSetMaxGPUMem( vcsbeam_context *vm, uintptr_t max_gpu_mem_bytes )
{
    vm->max_gpu_mem_bytes = max_gpu_mem_bytes;

    // Requested maximum can't be more that available memory
    struct cudaDeviceProp gpu_properties;
    cudaGetDeviceProperties( &gpu_properties, 0 );

    if (max_gpu_mem_bytes == 0) // Default behaviour: "0" = just use maximum available
    {
        vm->max_gpu_mem_bytes = gpu_properties.totalGlobalMem;
    }
    else if (max_gpu_mem_bytes > gpu_properties.totalGlobalMem )
    {
        fprintf( stderr, "warning: vmSetMaxGPUMem(): requested maximum (%lu) "
                "exceeds available memory (%lu). Setting to max available.\n",
                vm->max_gpu_mem_bytes, gpu_properties.totalGlobalMem );
        vm->max_gpu_mem_bytes = gpu_properties.totalGlobalMem;
    }

    // (This currently only accounts for some of the memory needed)
    vm->chunks_per_second = (vm->bytes_per_second + vm->Jv_size_bytes) /
        vm->max_gpu_mem_bytes + 1;

    // Make sure the number of chunks is divisible by the number of samples (per second)
    while ( vm->sample_rate % vm->chunks_per_second != 0 )
        vm->chunks_per_second++;

    // Calculate the amount of gpu memory needed
    vm->d_v_size_bytes = vm->bytes_per_second / vm->chunks_per_second;
    vm->d_Jv_size_bytes = vm->Jv_size_bytes / vm->chunks_per_second;
}

void vmPushChunk( vcsbeam_context *vm )
// Loads a "chunk" of data onto the GPU
{
    logger_start_stopwatch( vm->log, "upload", false );

    int chunk = vm->chunk_to_load % vm->chunks_per_second;
    char *ptrHost = (char *)vm->v->buffer + chunk*vm->d_v_size_bytes;

    cudaMemcpy( vm->d_v, ptrHost, vm->d_v_size_bytes, cudaMemcpyHostToDevice );
    cudaCheckErrors( "vmMemcpyNextChunk: cudaMemcpy failed" );

    logger_stop_stopwatch( vm->log, "upload" );
}


vm_error vmReadNextSecond( vcsbeam_context *vm )
{
    // Shorthand variables
    uintptr_t timestep_idx = vm->current_gps_idx;
    uintptr_t ntimesteps   = vm->num_gps_seconds_to_process;
    uint64_t  gps_second   = vm->gps_seconds_to_process[timestep_idx];

    // Make sure the buffer is allocated
    if (vm->v == NULL)
        return VM_READ_BUFFER_NOT_SET;

    // Make sure the buffer is not locked
    if (vm->v->locked)
        return VM_READ_BUFFER_LOCKED;

    // Return "end of data" if the timestep idx is too high
    if (timestep_idx >= ntimesteps)
        return VM_END_OF_DATA;

    // Now lock the buffer!
    vm->v->locked = true;

    int coarse_chan_idx = vm->coarse_chan_idxs_to_process[0];
    sprintf( vm->log_message, "--- Processing GPS second %ld [%lu/%lu], Coarse channel %lu [%d/%d] ---",
                gps_second,
                timestep_idx+1,
                ntimesteps,
                vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number,
                coarse_chan_idx,
                vm->num_coarse_chans_to_process );
    logger_timed_message( vm->log, vm->log_message );

    logger_start_stopwatch( vm->log, "read", true );

    vmReadBufferCopyMargin( vm->v );

    if (mwalib_voltage_context_read_second(
                vm->vcs_context,
                gps_second,
                1,
                coarse_chan_idx,
                vm->v->read_ptr,
                vm->v->read_size,
                vm->error_message,
                ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
    {
        fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s", vm->error_message );
        exit(EXIT_FAILURE);
    }

    logger_stop_stopwatch( vm->log, "read" );

    // Increment the count of number of seconds read
    vm->current_gps_idx++;

    return VM_SUCCESS;
}

void vmPushJ( vcsbeam_context *vm )
{
    cudaMemcpy( vm->d_J, vm->J, vm->J_size_bytes, cudaMemcpyHostToDevice );
    cudaCheckErrors( "vmMemcpyJ: cudaMemcpy failed" );
}


void vmCreateCudaStreams( vcsbeam_context *vm )
{
    vm->streams = (cudaStream_t *)malloc( vm->npointing * sizeof(cudaStream_t) );

    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        cudaStreamCreate( &(vm->streams[p]) );
        cudaCheckErrors( "vmCreateCudaStreams: cudaStreamCreate failed" );
    }
}

void vmDestroyCudaStreams( vcsbeam_context *vm )
{
    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        cudaStreamDestroy( vm->streams[p] );
        cudaCheckErrors( "vmDestroyCudaStreams: cudaStreamDestroy failed" );
    }

    free( vm->streams );
}

void vmCreateStatistics( vcsbeam_context *vm, mpi_psrfits *mpfs )
{
    uintptr_t nchan  = vm->nfine_chan;

    vm->offsets_size = vm->npointing*nchan*NSTOKES*sizeof(float);
    vm->scales_size  = vm->npointing*nchan*NSTOKES*sizeof(float);
    vm->Cscaled_size = vm->npointing*mpfs[0].coarse_chan_pf.sub.bytes_per_subint;

    cudaMalloc( (void **)&vm->d_offsets, vm->offsets_size );
    cudaCheckErrors( "vmCreateStatistics: cudaMalloc(offsets) failed" );
    cudaMalloc( (void **)&vm->d_scales,  vm->scales_size );
    cudaCheckErrors( "vmCreateStatistics: cudaMalloc(scales) failed" );
    cudaMalloc( (void **)&vm->d_Cscaled, vm->Cscaled_size );
    cudaCheckErrors( "vmCreateStatistics: cudaMalloc(Cscaled) failed" );

    cudaMallocHost( (void **)&vm->offsets, vm->offsets_size );
    cudaCheckErrors( "vmCreateStatistics: cudaMallocHost(offsets) failed" );
    cudaMallocHost( (void **)&vm->scales,  vm->scales_size );
    cudaCheckErrors( "vmCreateStatistics: cudaMallocHost(scales) failed" );
    cudaMallocHost( (void **)&vm->Cscaled, vm->Cscaled_size );
    cudaCheckErrors( "vmCreateStatistics: cudaMallocHost(Cscaled) failed" );
}

void vmDestroyStatistics( vcsbeam_context *vm )
{
    cudaFreeHost( vm->offsets );
    cudaCheckErrors( "vmDestroyStatistics: cudaFreeHost(offsets) failed" );
    cudaFreeHost( vm->scales );
    cudaCheckErrors( "vmDestroyStatistics: cudaFreeHost(scales) failed" );
    cudaFreeHost( vm->Cscaled );
    cudaCheckErrors( "vmDestroyStatistics: cudaFreeHost(Cscaled) failed" );

    cudaFree( vm->d_offsets );
    cudaCheckErrors( "vmDestroyStatistics: cudaFree(offsets) failed" );
    cudaFree( vm->d_scales );
    cudaCheckErrors( "vmDestroyStatistics: cudaFree(scales) failed" );
    cudaFree( vm->d_Cscaled );
    cudaCheckErrors( "vmDestroyStatistics: cudaFree(Cscaled) failed" );
}

void vmSetNumPointings( vcsbeam_context *vm, unsigned int npointings )
{
    uintptr_t npol      = vm->obs_metadata->num_ant_pols; // = 2
    vm->npointing       = npointings;
    vm->S_size_bytes    = vm->npointing * vm->nchan * NSTOKES * vm->sample_rate * sizeof(float);
    vm->d_S_size_bytes  = vm->S_size_bytes;
    vm->e_size_bytes    = vm->npointing * vm->sample_rate * vm->nchan * npol * sizeof(cuDoubleComplex);
    vm->d_e_size_bytes  = vm->e_size_bytes;
}

void vmCreateFilenames( vcsbeam_context *vm )
/* Create an array of filenames; free with vmDestroyFilenames()
 */
{
    // Shorthand variables
    unsigned long int begin_gps       = vm->gps_seconds_to_process[0];
    int               nseconds        = vm->num_gps_seconds_to_process;
    uintptr_t         coarse_chan_idx = vm->coarse_chan_idxs_to_process[0];
    int               ncoarse_chans   = vm->num_coarse_chans_to_process;

    // Check that the number of seconds is non-negative
    if (nseconds <= 0)
    {
        fprintf( stderr, "error: vmCreateFilenames: nseconds (= %d) must be >= 1\n",
                 nseconds );
        exit(EXIT_FAILURE);
    }

    // Check that the number of requested course channels does not exceed the number of available coarse channels
    if (coarse_chan_idx + ncoarse_chans > vm->obs_metadata->num_metafits_coarse_chans)
    {
        fprintf( stderr, "error: vmCreateFilenames: requested coarse channels higher than "
                "what is available\n" );
        exit(EXIT_FAILURE);
    }

    // Work out the number of files
    uint64_t end_gps = begin_gps + nseconds - 1;
    uint64_t file_start_timestep_second = begin_gps - (begin_gps % vm->seconds_per_file);
    uint64_t file_end_timestep_second   = end_gps - (end_gps % vm->seconds_per_file);
    vm->nfiletimes = (file_end_timestep_second - file_start_timestep_second)/vm->seconds_per_file + 1;
    vm->nfiles = ncoarse_chans * vm->nfiletimes;

    uint64_t gps_second;

    // Allocate memory for the file name list
    char filename[MAX_COMMAND_LENGTH]; // Just the mwalib-generated filename (without the path)
    vm->filenames = (char **)malloc( vm->nfiles * ncoarse_chans * sizeof(char *) ); // The full array of filenames, including the paths

    // Allocate memory and write filenames
    unsigned int t_idx, f_idx;
    uintptr_t c_idx;
    for (t_idx = 0; t_idx < vm->nfiletimes; t_idx++)
    {
        // Get the GPS second for this time index
        gps_second = file_start_timestep_second + t_idx*vm->seconds_per_file;

        for (c_idx = coarse_chan_idx; c_idx < coarse_chan_idx + ncoarse_chans; c_idx++)
        {
            //f_idx = second*ncoarse_chans + c_idx - coarse_chan_idx; // <-- this order seems to break mwalib (v0.9.4)
            f_idx = (c_idx - coarse_chan_idx)*vm->nfiletimes + t_idx;
            vm->filenames[f_idx] = (char *)malloc( MAX_COMMAND_LENGTH );
            memset( filename, 0, MAX_COMMAND_LENGTH );

            vmGetVoltFilename( vm, c_idx, gps_second, filename );

            sprintf( vm->filenames[f_idx], "%s/%s",
                    vm->datadir, filename );
        }
    }
}

/**
 * vmGetVoltFilename
 * =================
 *
 * Returns the filename for the given channel index and GPS second,
 * assuming the observation specified by VM->OBS_CONTEXT.
 * Expects VM->SECONDS_PER_FILE to be set (via vmLoadObsMetafits).
 * Result put into FILENAME (assumed already allocated).
 */
void vmGetVoltFilename( vcsbeam_context *vm, unsigned int coarse_chan_idx, uint64_t gps_second, char *filename )
{
    uint64_t t0_gps_second = vm->obs_metadata->metafits_timesteps[0].gps_time_ms/1000;
    uintptr_t timestep_idx = (gps_second - t0_gps_second) / vm->seconds_per_file;

    if (mwalib_metafits_get_expected_volt_filename(
                vm->obs_context,
                timestep_idx,
                coarse_chan_idx,
                filename,
                MAX_COMMAND_LENGTH,
                vm->error_message,
                ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): vmGetVoltFilename: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }
}

/**
 * vmGetLegacyVoltFilename
 * =======================
 *
 * Same as vmGetVoltFilename, but return the filename as if the observation
 * were really a legacy VCS observation.
 */
void vmGetLegacyVoltFilename( vcsbeam_context *vm, unsigned int coarse_chan_idx, uint64_t gps_second, char *filename )
{
    uint64_t t0_gps_second = vm->obs_metadata->metafits_timesteps[0].gps_time_ms/1000;
    uintptr_t timestep_idx = gps_second - t0_gps_second;

    if (mwalib_metafits_get_expected_volt_filename(
                vm->obs_context_legacy,
                timestep_idx,
                coarse_chan_idx,
                filename,
                MAX_COMMAND_LENGTH,
                vm->error_message,
                ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): vmGetVoltFilename: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }
}

void vmDestroyFilenames( vcsbeam_context *vm )
{
    if (vm->filenames != NULL)
    {
        int i;
        for (i = 0; i < vm->nfiles; i++)
            free( vm->filenames[i] );
        free( vm->filenames );
    }

    vm->nfiles = 0;
}


void vmLoadObsMetafits( vcsbeam_context *vm, char *filename )
{
    // Create OBS_CONTEXT
    if (mwalib_metafits_context_new2( filename, &vm->obs_context, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_CONTEXT_LEGACY
    if (mwalib_metafits_context_new( filename, VCSLegacyRecombined, &vm->obs_context_legacy,
            vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create legacy metafits context: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( vm->obs_context, NULL, NULL, &vm->obs_metadata, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    vm->seconds_per_file = (vm->obs_metadata->mwa_version == VCSMWAXv2 ? 8 : 1);
                        // ^-- Raise an issue with mwalib (I shouldn't have to do this)
}

void vmLoadCalMetafits( vcsbeam_context *vm, char *filename )
{
    // Create OBS_CONTEXT
    if (mwalib_metafits_context_new2( filename, &vm->cal_context, vm->error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( vm->cal_context, NULL, NULL, &vm->cal_metadata, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }
}

void vmGetVoltageMetadata( vcsbeam_context *vm )
/* Create the voltage metadata structs using mwalib's API.
 */
{
    // Shorthand variables
    unsigned long int begin_gps       = vm->gps_seconds_to_process[0];
    int               nseconds        = vm->num_gps_seconds_to_process;
    int               ncoarse_chans   = vm->num_coarse_chans_to_process;

    // If nseconds is an invalid number (<= 0), make it equal to the max possible for this obs
    if (nseconds <= 0)
    {
        nseconds = vm->obs_metadata->num_metafits_timesteps - (begin_gps - vm->obs_metadata->metafits_timesteps[0].gps_time_ms/1000);
    }

    // Ditto for the ncoarse_chans
    if (ncoarse_chans < 1)
    {
        fprintf( stderr, "error: vmGetVoltageMetadata: number of coarse chans must be >1\n" );
        exit(EXIT_FAILURE);
    }

    // Create list of filenames
    vmCreateFilenames( vm );

    // Create an mwalib voltage context, voltage metadata, and new obs metadata (now with correct antenna ordering)
    // (MWALIB is expecting a const array, so we will give it one!)
    const char **voltage_files = (const char **)malloc( sizeof(char *) * vm->nfiles );
    int i;
    for (i = 0; i < vm->nfiles; i++)
    {
        voltage_files[i] = vm->filenames[i];
    }

    // Create VCS_CONTEXT
    if (mwalib_voltage_context_new( vm->obs_metadata->metafits_filename, voltage_files, vm->nfiles, &vm->vcs_context, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create voltage context: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Create VCS_METADATA
    if (mwalib_voltage_metadata_get( vm->vcs_context, &vm->vcs_metadata, vm->error_message, ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot get metadata: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Replace the existing OBS_METADATA with a new one that knows this is a VCS observation
    mwalib_metafits_metadata_free( vm->obs_metadata );
    if (mwalib_metafits_metadata_get( NULL, NULL, vm->vcs_context, &vm->obs_metadata, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata from voltage context: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Free memory
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


void vmSetNumNotFlaggedRFInputs( vcsbeam_context *vm )
{
    vm->num_not_flagged = 0;
    uintptr_t i;
    for (i = 0; i < vm->obs_metadata->num_rf_inputs; i++)
        if (!(vm->obs_metadata->rf_inputs[i].flagged))
            vm->num_not_flagged++;
}

Rfinput *find_matching_rf_input( MetafitsMetadata *metadata, Rfinput *rfinput )
{
    // Find the input in METADATA that has the matching tile_id and polarisation
    uintptr_t i;
    char pol;
    uint32_t tile_id;
    for (i = 0; i < metadata->num_rf_inputs; i++)
    {
        tile_id   = metadata->rf_inputs[i].tile_id;
        pol       = *(metadata->rf_inputs[i].pol);

        if (rfinput->tile_id == tile_id && *(rfinput->pol) == pol)
        {
            return &(metadata->rf_inputs[i]);
        }
    }

    return NULL;
}


Antenna *find_matching_antenna( MetafitsMetadata *metadata, Rfinput *rfinput )
{
    // Find the input in METADATA that has the matching tile_id and polarisation
    uintptr_t a;
    uint32_t tile_id;
    for (a = 0; a < metadata->num_ants; a++)
    {
        tile_id   = metadata->antennas[a].tile_id;

        if (rfinput->tile_id == tile_id)
        {
            return &(metadata->antennas[a]);
        }
    }

    return NULL;
}

void get_mwalib_version( char *version_str )
/* Assumes that version_str is already allocated, and is big enough
 */
{
    sprintf( version_str, "%u.%u.%u",
            mwalib_get_version_major(),
            mwalib_get_version_minor(),
            mwalib_get_version_patch()
           );
}


void vmParsePointingFile( vcsbeam_context *vm, const char *filename )
/* Parse the given file in FILENAME and create arrays of RAs and DECs.
 * This function allocates memory for ras_hours and decs_degs arrays.
 * Will be destroyed when vcsbeam context is destroyed.
 */
{
    // Print a log message
    sprintf( vm->log_message, "Reading pointings file %s", filename );
    logger_timed_message( vm->log, vm->log_message );

    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL)
    {
        fprintf( stderr, "error: cannot open pointings file %s\n", filename );
        exit(EXIT_FAILURE);
    }

    // Do one pass through the file to count "words"
    // The RAs and Decs are expected to be whitespace-delimited
    int nwords = 0;
    char word[64];
    while (fscanf( f, "%s", word ) != EOF)
        nwords++;

    // Check that we have an even number of words (they should be in RA/Dec pairs)
    if (nwords % 2 != 0)
    {
        fprintf( stderr, "error: cannot parse pointings file %s\n", filename );
        exit(EXIT_FAILURE);
    }
    vmSetNumPointings( vm, nwords/2 );

    // Allocate memory
    vm->ras_hours = (double *)malloc( vm->npointing * sizeof(double) );
    vm->decs_degs = (double *)malloc( vm->npointing * sizeof(double) );

    // Rewind to beginning of file, read the words in again, and parse them
    rewind( f );
    char ra_str[64], dec_str[64];
    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        // Read in the next Ra/Dec pair
        fscanf( f, "%s %s", ra_str, dec_str );

        // Parse them and make them decimal
        vm->ras_hours[p] = parse_ra( ra_str );
        vm->decs_degs[p] = parse_dec( dec_str );
    }

    // Close the file
    fclose( f );

    // Set up CUDA streams (one stream per pointing)
    vmCreateCudaStreams( vm );
}

void vmReportPerformanceStats( vcsbeam_context *vm )
{
    logger_report_all_stats( vm->log );
}

void vmPrintTitle( vcsbeam_context *vm, const char *title )
{
    sprintf( vm->log_message, "------- VCSBeam (%s): %s -------",
            VCSBEAM_VERSION, title );
    logger_message( vm->log, vm->log_message );
}
