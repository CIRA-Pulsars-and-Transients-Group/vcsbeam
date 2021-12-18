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
        vm->mpi_rank = 0;
    }
    vm->writer = 0;

    // TODO: Change this to give user flexibility of how to use mpi structure
    vm->ncoarse_chans   = vm->mpi_size;
    vm->coarse_chan_idx = vm->mpi_rank;

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

    vm->v_size_bytes          = 0;
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

    // Default: data is legacy VCS format (VB_INT4)
    vm->datatype = VB_INT4;

    // Calibration
    init_calibration( &vm->cal );

    // Initially not bound to any observations
    vm->obs_context  = NULL;
    vm->obs_metadata = NULL;

    vm->vcs_context  = NULL;
    vm->vcs_metadata = NULL;

    vm->cal_context  = NULL;
    vm->cal_metadata = NULL;

    // Start a logger
    vm->log = create_logger( stdout, vm->mpi_rank );
    logger_add_stopwatch( vm->log, "read", "Reading in data" );
    logger_add_stopwatch( vm->log, "delay", "Calculating geometric and cable delays" );
    logger_add_stopwatch( vm->log, "calc", "Calculating tied-array beam" );
    logger_add_stopwatch( vm->log, "ipfb", "Inverting the PFB" );
    logger_add_stopwatch( vm->log, "splice", "Splicing coarse channels together" );
    logger_add_stopwatch( vm->log, "write", "Writing out data to file" );

    // Initialise pointing RAs and Decs to NULL
    vm->ras_hours = NULL;
    vm->decs_degs = NULL;

    // Return the new struct pointer
    return vm;
}

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
            fprintf( stderr, "vmBindObsData: error: "
                    "this observation does not appear to be a VCS observation\n" );
            exit(EXIT_FAILURE);
            break;
    }

    // Compute the "shorthand" variables
    vm->sample_rate = vm->vcs_metadata->num_samples_per_voltage_block *
                      vm->vcs_metadata->num_voltage_blocks_per_second;

    vm->bytes_per_second = vm->vcs_metadata->num_voltage_blocks_per_second *
                           vm->vcs_metadata->voltage_block_size_bytes;

    vm->nchan = vm->obs_metadata->num_volt_fine_chans_per_coarse;

    // Assume that the whole GPU is available
    uintptr_t nant      = vm->obs_metadata->num_ants;
    uintptr_t nvispol   = vm->obs_metadata->num_visibility_pols; // = 4 (PP, PQ, QP, QQ)
    uintptr_t npol      = vm->obs_metadata->num_ant_pols; // = 2

    vm->v_size_bytes    = vm->bytes_per_second;
    vm->d_v_size_bytes  = vm->v_size_bytes;
    vm->pol_idxs_size_bytes   = nant * sizeof(uint32_t);
    vm->d_pol_idxs_size_bytes = vm->pol_idxs_size_bytes;
    vm->D_size_bytes    = nant * vm->nchan * nvispol * sizeof(cuDoubleComplex);
    vm->d_D_size_bytes  = vm->D_size_bytes;
    vm->J_size_bytes    = nant * vm->nchan * nvispol * sizeof(cuDoubleComplex);
    vm->d_J_size_bytes  = vm->J_size_bytes;
    vm->e_size_bytes    = vm->npointing * vm->sample_rate * vm->nchan * npol * sizeof(cuDoubleComplex);
    vm->d_e_size_bytes  = vm->e_size_bytes;
    vm->S_size_bytes    = vm->npointing * vm->nchan * NSTOKES * vm->sample_rate * sizeof(float);
    vm->d_S_size_bytes  = vm->S_size_bytes;
    vm->Jv_size_bytes   = nant * vm->nchan * vm->sample_rate * sizeof(cuDoubleComplex);
    vm->d_Jv_size_bytes = vm->Jv_size_bytes;

}


void vmBindCalData( vcsbeam_context *vm,
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

    if (vm->cal.cal_type == CAL_RTS)
    {
        vmLoadRTSSolution( vm, vm->cal.use_bandpass, vm->cal.caldir, vm->coarse_chan_idx );
    }
    else if (vm->cal.cal_type == CAL_OFFRINGA)
    {
        vmLoadOffringaSolution( vm, vm->coarse_chan_idx, vm->cal.caldir );
    }

    // Flag antennas that need flagging
    vmSetCustomTileFlags( vm, flags_file, &vm->cal );
}

Antenna *find_antenna_by_name( MetafitsMetadata *obs_metadata, char *tile_name )
/* Returns the index (to obs_metadata->antennas[]) for the antenna with the
 * specified name. If no antenna with that name exists, returns
 * NO_ANTENNA_FOUND (= -1)
 */
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


void destroy_vcsbeam_context( vcsbeam_context *vm )
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

    destroy_logger( vm->log );

    // Free the RA and Dec pointing arrays
    if (vm->ras_hours != NULL)  free( vm->ras_hours );
    if (vm->decs_degs != NULL)  free( vm->decs_degs );

    // Calibration
    free_calibration( &vm->cal );

    // Finalise MPI
    if (vm->use_mpi)
        MPI_Finalize();

    // Finally, free the struct itself
    free( vm );
}

void set_vcsbeam_fine_output( vcsbeam_context *vm, bool switch_on )
/* Turns on/off fine channelised output
 */
{
    vm->output_fine_channels = switch_on;
    if (vm->obs_metadata->mwa_version == VCSMWAXv2)
        vm->do_forward_pfb = switch_on;
}

void set_vcsbeam_coarse_output( vcsbeam_context *vm, bool switch_on )
/* Turns on/off coarse channelised output
 */
{
    vm->output_coarse_channels = switch_on;
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
        vm->do_inverse_pfb = switch_on;
}

void vmMallocVHost( vcsbeam_context *vm )
{
    cudaMallocHost( (void **)&(vm->v), vm->v_size_bytes );
    cudaCheckErrors( "vmMallocVHost: cudaMallocHost failed" );
}

void vmFreeVHost( vcsbeam_context *vm )
{
    cudaFreeHost( vm->v );
    cudaCheckErrors( "vmFreeVHost: cudaFreeHost failed" );
}

void vmMallocVDevice( vcsbeam_context *vm )
{
    cudaMalloc( (void **)&vm->d_v,  vm->d_v_size_bytes );
    cudaCheckErrors( "vmMallocVDevice: cudaMalloc failed" );
}

void vmFreeVDevice( vcsbeam_context *vm )
{
    cudaFree( vm->d_v );
    cudaCheckErrors( "vmFreeVDevice: cudaFree failed" );
}

void vmMallocJVHost( vcsbeam_context *vm )
{
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
    vm->chunks_per_second = (vm->v_size_bytes + vm->Jv_size_bytes) /
        vm->max_gpu_mem_bytes + 1;

    // Make sure the number of chunks is divisible by the number of samples (per second)
    while ( vm->sample_rate % vm->chunks_per_second != 0 )
        vm->chunks_per_second++;

    // Calculate the amount of gpu memory needed
    vm->d_v_size_bytes = vm->v_size_bytes / vm->chunks_per_second;
    vm->d_Jv_size_bytes = vm->Jv_size_bytes / vm->chunks_per_second;
}

void vmPushChunk( vcsbeam_context *vm )
{
    // Loads a "chunk" of data onto the GPU
    int chunk = vm->chunk_to_load % vm->chunks_per_second;
    char *ptrHost = (char *)vm->v + chunk*vm->d_v_size_bytes;
    cudaMemcpy( vm->d_v, ptrHost, vm->d_v_size_bytes, cudaMemcpyHostToDevice );
    cudaCheckErrors( "vmMemcpyNextChunk: cudaMemcpy failed" );
}


vcsbeam_error vmReadNextSecond( vcsbeam_context *vm )
{
    uintptr_t timestep_idx = vm->chunk_to_load / vm->chunks_per_second;
    uintptr_t ntimesteps   = vm->num_gps_seconds_to_process;
    uint64_t  gps_second   = vm->gps_seconds_to_process[timestep_idx];

    // Return "end of data" if the timestep idx is too high
    if (timestep_idx >= ntimesteps)
        return VB_EOD;

    sprintf( vm->log_message, "--- Processing GPS second %ld [%lu/%lu], Coarse channel %lu [%d/%d] ---",
                gps_second, timestep_idx+1, ntimesteps,
                vm->obs_metadata->metafits_coarse_chans[vm->coarse_chan_idxs_to_process[0]].rec_chan_number,
                vm->coarse_chan_idxs_to_process[0], vm->num_coarse_chans_to_process );
    logger_message( vm->log, vm->log_message );

    logger_start_stopwatch( vm->log, "read", true );

    if (mwalib_voltage_context_read_second(
                vm->vcs_context,
                vm->gps_seconds_to_process[timestep_idx],
                1,
                vm->coarse_chan_idxs_to_process[0],
                vm->v,
                vm->v_size_bytes,
                vm->error_message,
                ERROR_MESSAGE_LEN ) != EXIT_SUCCESS)
    {
        fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s", vm->error_message );
        exit(EXIT_FAILURE);
    }

    logger_stop_stopwatch( vm->log, "read" );

    return VB_SUCCESS;
}

void vmPushJ( vcsbeam_context *vm )
{
    cudaMemcpy( vm->d_J, vm->J, vm->J_size_bytes, cudaMemcpyHostToDevice );
    cudaCheckErrors( "vmMemcpyJ: cudaMemcpy failed" );
}


void vmCreatePrimaryBeam( vcsbeam_context *vm )
{
    create_primary_beam( &vm->pb, vm->obs_metadata, vm->coarse_chan_idxs_to_process[0], vm->npointing );
}

void vmCreateGeometricDelays( vcsbeam_context *vm )
{
    create_geometric_delays( &vm->gdelays, vm->obs_metadata, vm->vcs_metadata, vm->coarse_chan_idxs_to_process[0], vm->npointing );
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
    uintptr_t nchan  = vm->obs_metadata->num_volt_fine_chans_per_coarse;

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


void vmLoadObsMetafits( vcsbeam_context *vm, char *filename )
{
    // Create OBS_CONTEXT
    if (mwalib_metafits_context_new2( filename, &vm->obs_context, vm->error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Create OBS_METADATA
    if (mwalib_metafits_metadata_get( vm->obs_context, NULL, NULL, &vm->obs_metadata, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }
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

