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

#include "gpu_macros.h"

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
 * @param vm The VCSBeam context struct
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
    vm->cal_coarse_chan_idxs_to_process = (int *)malloc( vm->num_coarse_chans_to_process * sizeof(int) );

    int g;
    for (g = 0; g < vm->num_gps_seconds_to_process; g++)
    {
        vm->gps_seconds_to_process[g] = g + first_gps_second;
    }

    int c;
    for (c = 0; c < num_coarse_chans_to_process; c++)
    {
        vm->coarse_chan_idxs_to_process[c] = c + first_coarse_chan_idx;
        vm->cal_coarse_chan_idxs_to_process[c] = c + first_coarse_chan_idx;
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
 * @param vm The VCSBeam context struct
 * @param caldir The directory containing RTS solution files, OR the path of
 *        an Offringa-style calibration solution file
 * @param cal_type Either `CAL_RTS` or `CAL_OFFRINGA`
 * @param use_bandpass Whether to include the Bandpass information (relevant
 *        for RTS solutions only)
 * @param flags_file A file containing names of (extra) tiles to be flagged,
 *        or `NULL`
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

    if (flags_file != NULL)
    {
        vm->cal.flags_file   = (char *)malloc( strlen( flags_file ) + 1 );
        strcpy( vm->cal.flags_file, flags_file );
    }
    else
        vm->cal.flags_file = NULL;
}

/**
 * Reads in a calibration solution.
 *
 * @param vm The VCSBeam context struct
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
 * @param vm The VCSBeam context struct
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
    free( vm->cal_coarse_chan_idxs_to_process );

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
 * @param vm The VCSBeam context struct
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

/**
 * Allocates memory for the input voltages on the CPU.
 *
 * @param vm The VCSBeam context struct
 *
 * The memory is allocated as a `host_buffer` struct (see vmInitReadBuffer()
 * for a full description).
 * The amount of memory allocated depends on whether the observation is
 * Legacy or MWAX.
 *
 * If the input is legacy, the read buffer can be just enough
 * for one second's worth of data.
 *
 * However, for MWAX, the read buffer must be slightly bigger
 * to accommodate the PFB's extra taps.
 * (For now, this is fixed to one voltage block. Changing this
 * will break the gpu kernels that do the fine PFB, for which
 * the offset into the data is currently hard-coded.)
 */
void vmMallocVHost( vcsbeam_context *vm )
{
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
        vm->v = vmInitReadBuffer( vm->bytes_per_second, 0 );
    else if (vm->obs_metadata->mwa_version == VCSMWAXv2)
        vm->v = vmInitReadBuffer( vm->bytes_per_second, vm->vcs_metadata->voltage_block_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf J}{\bf v}\f$ on the CPU.
 *
 * @param vm The VCSBeam context struct
 *
 * Pointers to the newly allocated memory are given in
 * `vm&rarr;Jv_P` and `vm&rarr;Jv_Q`.
 */
void vmMallocJVHost( vcsbeam_context *vm )
{
    vm->Jv_size_bytes =
        vm->npointing *
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->fine_sample_rate *
        sizeof(gpuDoubleComplex);

    gpuMallocHost( (void **)&(vm->Jv_P), vm->Jv_size_bytes );
    gpuMallocHost( (void **)&(vm->Jv_Q), vm->Jv_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf e}\f$ on the CPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;e`.
 */
void vmMallocEHost( vcsbeam_context *vm )
{
    vm->e_size_bytes =
        vm->npointing *
        vm->fine_sample_rate *
        vm->nfine_chan *
        vm->obs_metadata->num_ant_pols *
        sizeof(gpuDoubleComplex);

    // Allocate memory on host
    gpuMallocHost( (void **)&(vm->e), vm->e_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf S}\f$ on the CPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;S`.
 */
void vmMallocSHost( vcsbeam_context *vm )
{
    vm->S_size_bytes = vm->npointing * vm->nfine_chan * vm->out_nstokes * vm->fine_sample_rate * sizeof(float);

    // Allocate memory on host
    gpuMallocHost( (void **)&(vm->S), vm->S_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf J}\f$ on the CPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;J`.
 */
void vmMallocJHost( vcsbeam_context *vm )
{
    vm->J_size_bytes =
        vm->npointing *
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(gpuDoubleComplex);

    // Allocate memory on device
    gpuMallocHost( (void **)&(vm->J), vm->J_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf D}\f$ on the CPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;D`.
 */
void vmMallocDHost( vcsbeam_context *vm )
{
    vm->D_size_bytes =
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(gpuDoubleComplex);

    // Allocate memory on device
    gpuMallocHost( (void **)&(vm->D), vm->D_size_bytes );
}

/**
 * Allocates memory for the polarisation indexes on the CPU.
 *
 * @param vm The VCSBeam context struct
 *
 * Pointers to the newly allocated memory are given in
 * `vm&rarr;polP_idxs` and `vm&rarr;polQ_idxs`.
 */
void vmMallocPQIdxsHost( vcsbeam_context *vm )
{
    vm->pol_idxs_size_bytes = vm->obs_metadata->num_ants * sizeof(uint32_t);

    // Allocate memory on device
    gpuMallocHost( (void **)&(vm->polP_idxs), vm->pol_idxs_size_bytes );
    gpuMallocHost( (void **)&(vm->polQ_idxs), vm->pol_idxs_size_bytes );
}

/**
 * Frees memory allocated for the input voltages on the CPU.
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeVHost( vcsbeam_context *vm )
{
    vmFreeReadBuffer( vm->v );
}

/**
 * Frees memory allocated with vmMallocJVHost().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeJVHost( vcsbeam_context *vm )
{
    gpuHostFree( vm->Jv_P );
    gpuHostFree( vm->Jv_Q );
}

/**
 * Frees memory allocated with vmMallocEHost().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeEHost( vcsbeam_context *vm )
{
    gpuHostFree( vm->e );
}

/**
 * Frees memory allocated with vmMallocSHost().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeSHost( vcsbeam_context *vm )
{
    gpuHostFree( vm->S );
}

/**
 * Frees memory allocated with vmMallocJHost().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeJHost( vcsbeam_context *vm )
{
    gpuHostFree( vm->J );
}

/**
 * Frees memory allocated with vmMallocDHost().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeDHost( vcsbeam_context *vm )
{
    gpuHostFree( vm->D );
}

/**
 * Frees memory allocated with vmMallocPQIdxsHost().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreePQIdxsHost( vcsbeam_context *vm )
{
    gpuHostFree( vm->polP_idxs );
    gpuHostFree( vm->polQ_idxs );
}

/**
 * Allocates memory for the input voltages on the GPU.
 *
 * @param vm The VCSBeam context struct
 *
 * If the observation is Legacy, then this function allocates new memory and
 * sets `vm&rarr;d_v_size_bytes` to the address of the new memory block.
 *
 * If the observation is MWAX, then this function only copies the value of
 * `vm&rarr;fpfb&rarr;d_vcs_data` to `vm&rarr;d_v_size_bytes`, (without
 * checking whether it points to valid GPU memory).
 *
 * Only one "chunk" of memory is allocated (where each second of input data is
 * divided up into one or more chunks, depending on the amount of memory
 * available on the GPU).
 */
void vmMallocVDevice( vcsbeam_context *vm )
{
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
    {
        vm->d_v_size_bytes = vm->bytes_per_second / vm->chunks_per_second;
        gpuMalloc( (void **)&vm->d_v,  vm->d_v_size_bytes );
    }
    else // if (vm->obs_metadata->mwa_version == VCSMWAXv2)
    {
        vm->d_v            = vm->fpfb->d_vcs_data;
        vm->d_v_size_bytes = vm->fpfb->d_vcs_size;
        // TODO: ^^^ Check that these are already the same thing
    }
}

/**
 * Allocates memory for the quantities \f${\bf J}{\bf v}\f$ on the GPU.
 *
 * @param vm The VCSBeam context struct
 *
 * Pointers to the newly allocated memory are given in
 * `vm&rarr;d_Jv_P` and `vm&rarr;d_Jv_Q`.
 *
 * Only one "chunk" of memory is allocated (where each second of input data is
 * divided up into one or more chunks, depending on the amount of memory
 * available on the GPU).
 */
void vmMallocJVDevice( vcsbeam_context *vm )
{
    vm->d_Jv_size_bytes = 
        vm->npointing *
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->fine_sample_rate *
        sizeof(gpuDoubleComplex) /
        vm->chunks_per_second;

#ifdef DEBUG
    printf("%lu\n", vm->d_Jv_size_bytes);
    size_t mf, ma;
    cudaMemGetInfo(&mf, &ma); //TODO: This will break, need equiv. definitions
    printf("free: %zu ... total: %zu\n", mf, ma);
#endif

    gpuMalloc( (void **)&vm->d_Jv_P,  vm->d_Jv_size_bytes );

#ifdef DEBUG
    cudaMemGetInfo(&mf, &ma); //TODO: This will break, need equiv. definitions
    printf("free: %zu ... total: %zu\n", mf, ma);
#endif

    gpuMalloc( (void **)&vm->d_Jv_Q,  vm->d_Jv_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf e}\f$ on the GPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;d_e`.
 */
void vmMallocEDevice( vcsbeam_context *vm )
{
    vm->d_e_size_bytes =
        vm->npointing *
        vm->fine_sample_rate *
        vm->nfine_chan *
        vm->obs_metadata->num_ant_pols *
        sizeof(gpuDoubleComplex);

    // Allocate memory on device
    gpuMalloc( (void **)&(vm->d_e), vm->d_e_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf S}\f$ on the GPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;d_S`.
 */
void vmMallocSDevice( vcsbeam_context *vm )
{
    vm->d_S_size_bytes = vm->npointing * vm->nfine_chan * vm->out_nstokes * vm->fine_sample_rate * sizeof(float);

    // Allocate memory on device
    gpuMalloc( (void **)&(vm->d_S), vm->d_S_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf J}\f$ on the GPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;d_J`.
 */
void vmMallocJDevice( vcsbeam_context *vm )
{
    vm->d_J_size_bytes =
        vm->npointing *
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(gpuDoubleComplex);

    // Allocate memory on device
    gpuMalloc( (void **)&(vm->d_J), vm->d_J_size_bytes );
}

/**
 * Allocates memory for the quantity \f${\bf D}\f$ on the GPU.
 *
 * @param vm The VCSBeam context struct
 *
 * A pointer to the newly allocated memory is given in `vm&rarr;d_D`.
 */
void vmMallocDDevice( vcsbeam_context *vm )
{
    vm->d_D_size_bytes =
        vm->obs_metadata->num_ants *
        vm->nfine_chan *
        vm->obs_metadata->num_visibility_pols *
        sizeof(gpuDoubleComplex);

    // Allocate memory on device
    gpuMalloc( (void **)&(vm->d_D), vm->d_D_size_bytes );
}

/**
 * Allocates memory for the polarisation indexes on the GPU.
 *
 * @param vm The VCSBeam context struct
 *
 * Pointers to the newly allocated memory are given in
 * `vm&rarr;d_polP_idxs` and `vm&rarr;d_polQ_idxs`.
 */
void vmMallocPQIdxsDevice( vcsbeam_context *vm )
{
    vm->d_pol_idxs_size_bytes = vm->obs_metadata->num_ants * sizeof(uint32_t);

    // Allocate memory on device
    gpuMalloc( (void **)&(vm->d_polP_idxs), vm->d_pol_idxs_size_bytes );
    gpuMalloc( (void **)&(vm->d_polQ_idxs), vm->d_pol_idxs_size_bytes );
}

/**
 * Frees the GPU memory allocated with vmMallocVDevice().
 *
 * @param vm The VCSBeam context struct
 *
 * If the observation is Legacy, then free the memory.
 *
 * If the observation is MWAX, then do nothing; the memory should be freed via
 * a call to vmFreeForwardPFB().
 */
void vmFreeVDevice( vcsbeam_context *vm )
{
    if (vm->obs_metadata->mwa_version == VCSLegacyRecombined)
    {
        gpuFree( vm->d_v );
    }
}

/**
 * Frees memory allocated with vmMallocJVDevice().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeJVDevice( vcsbeam_context *vm )
{
    gpuFree( vm->d_Jv_P );
    gpuFree( vm->d_Jv_Q );
}

/**
 * Frees memory allocated with vmMallocEDevice().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeEDevice( vcsbeam_context *vm )
{
    gpuFree( vm->d_e );
}

/**
 * Frees memory allocated with vmMallocSDevice().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeSDevice( vcsbeam_context *vm )
{
    gpuFree( vm->d_S );
}

/**
 * Frees memory allocated with vmMallocJDevice().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeJDevice( vcsbeam_context *vm )
{
    gpuFree( vm->d_J );
}


/**
 * Frees memory allocated with vmMallocDDevice().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreeDDevice( vcsbeam_context *vm )
{
    gpuFree( vm->d_D );
}


/**
 * Frees memory allocated with vmMallocPQIdxsDevice().
 *
 * @param vm The VCSBeam context struct
 */
void vmFreePQIdxsDevice( vcsbeam_context *vm )
{
    gpuFree( vm->d_polP_idxs );
    gpuFree( vm->d_polQ_idxs );
}

/**
 * Tries to cleverly figure out how many chunks are needed to
 * fit everything on the GPU.
 *
 * @param vm The VCSBeam context struct
 * @param max_gpu_mem_bytes The maximum amount of GPU memory (in bytes) to use
 *
 * LOGIC IS CURRENTLY FAULTY AND INCOMPLETE. DO NOT USE!
 */
void vmSetMaxGPUMem( vcsbeam_context *vm, int nchunks )
{
    vm->chunks_per_second = nchunks;

    // Make sure the number of chunks is divisible by the number of samples (per second)
    if ( vm->sample_rate % vm->chunks_per_second != 0 )
    {
        logger_timed_message( vm->log, "Bad number of chunks requested: must divide number of samples per second of data." );
        exit(EXIT_FAILURE);
    }

    // Calculate the amount of gpu memory needed
    vm->d_v_size_bytes = vm->bytes_per_second / vm->chunks_per_second;
    vm->d_Jv_size_bytes = vm->Jv_size_bytes / vm->chunks_per_second;
}

/**
 * Loads a "chunk" of input data onto the GPU
 *
 * @param vm The VCSBeam context struct
 */
void vmPushChunk( vcsbeam_context *vm )
{
    logger_start_stopwatch( vm->log, "upload", false );

    int chunk = vm->chunk_to_load % vm->chunks_per_second;
    char *ptrHost = (char *)vm->v->buffer + chunk*vm->d_v_size_bytes;

    gpuMemcpy( vm->d_v, ptrHost, vm->d_v_size_bytes, gpuMemcpyHostToDevice );

    logger_stop_stopwatch( vm->log, "upload" );
}


/**
    @Cristian: here is a simple test to see if we can read as much data as we can in one go.
    And whether that will improve the I/O performance of the program. At first, it is all
    sequential. That is, a buffer is allocated, and populated when emptied, and consumed by
    vmReadNextSecond transparently.
*/

struct {
    char *data;
    uint64_t current_gps_second;
    uint64_t n_remaining_gps_seconds;
    uint64_t count;
    uint64_t size;
    uint64_t capacity;
} seconds_buffer = {NULL, 0u, 0u, 0u, 0u, 0u};


void read_from_buffer(vcsbeam_context *vm,
                    unsigned long gps_second_start,
                    size_t gps_second_count,
                    size_t voltage_coarse_chan_index,
                    signed char *buffer_ptr,
                    size_t buffer_len,
                    const char *error_message,
                    size_t error_message_length){
    
    if(seconds_buffer.data == NULL){
        // Initialise the structure
        //unsigned int nfiletimes;          // The number of "file" timesteps
        seconds_buffer.n_remaining_gps_seconds = vm->nfiletimes * vm->seconds_per_file;
        printf("CDP DEBUG: n_remaining_gps_seconds is set to: %lu, "
            "with nfiletimes = %lu and seconds per file = %lu\n",
            seconds_buffer.n_remaining_gps_seconds, vm->nfiletimes, vm->seconds_per_file);
        // the following must be greater than the number of seconds
        // in each file.
        // TODO: check gps_second_count is less than buffer size
        // TODO: find a clever way of doing this.
        size_t desired_seconds_in_buffer = 64;
        // this will ensure each file is read in full
        size_t total_seconds = vm->seconds_per_file * (desired_seconds_in_buffer / vm->seconds_per_file);
        size_t total_bytes = vm->bytes_per_second * total_seconds;
        seconds_buffer.capacity = total_seconds;
        printf("Will allocate %.4f GiB for 'seconds_buffer', corresponding to %lu seconds.\n", (total_bytes / (1024.0f * 1024.0f * 1024.0f)), total_seconds);
        seconds_buffer.data = (char*) malloc(total_bytes);
        seconds_buffer.current_gps_second = gps_second_start;
        if(!seconds_buffer.data){
            fprintf(stderr, "Error allocating memory for 'seconds_buffer'.\n");
            exit(1);
        }
    }
    // TODO: do not support gps_second_count != 1, because the mechanism to refill the buffer
    // might be more complicated
    if(seconds_buffer.count == seconds_buffer.size){
        // buffer is empty, refill it
        // TODO: check for read sizes less than buffer size
        // TODO: use read file function in mwalib next
        size_t seconds_to_read = seconds_buffer.capacity < seconds_buffer.n_remaining_gps_seconds ? \
            seconds_buffer.capacity : seconds_buffer.n_remaining_gps_seconds;
        size_t start_gps = seconds_buffer.current_gps_second;
        size_t end_gps = seconds_buffer.current_gps_second + seconds_to_read;
        size_t read_size = vm->bytes_per_second * vm->seconds_per_file;

        #pragma omp parallel for schedule(static) num_threads((seconds_to_read/vm->seconds_per_file))
        for(size_t sec_idx = start_gps; sec_idx < end_gps; sec_idx += vm->seconds_per_file){
            // We are assuming that the seconds_to_read is always a integer multiple of seconds_per_file!
            if(mwalib_voltage_context_read_second(vm->vcs_context, sec_idx, vm->seconds_per_file,
                    voltage_coarse_chan_index, seconds_buffer.data + read_size * (sec_idx - start_gps), read_size,
                    vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS){
                fprintf( stderr, "error: mwalib_voltage_context_read_file failed: %s", vm->error_message );
                exit(EXIT_FAILURE);
            }
        }
        seconds_buffer.count = 0u;
        seconds_buffer.size = seconds_to_read;
        seconds_buffer.n_remaining_gps_seconds -= seconds_to_read;
    }

    char *source = seconds_buffer.data + seconds_buffer.count * vm->bytes_per_second;
    memcpy(buffer_ptr, source, vm->bytes_per_second * gps_second_count);
    seconds_buffer.count += gps_second_count;
    // we are assuming seconds are read sequentially!
    seconds_buffer.current_gps_second += gps_second_count;
}


/**
 * Reads a second's worth of input data from the observation files.
 *
 * @param vm The VCSBeam context struct
 */
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

    /*
    mwalib_voltage_context_display(
            vm->vcs_context,
            vm->error_message,
            ERROR_MESSAGE_LEN);
    */
    read_from_buffer(
                vm,
                gps_second,
                1,
                coarse_chan_idx,
                vm->v->read_ptr,
                vm->v->read_size,
                vm->error_message,
                ERROR_MESSAGE_LEN);

    logger_stop_stopwatch( vm->log, "read" );

    // Increment the count of number of seconds read
    vm->current_gps_idx++;

    return VM_SUCCESS;
}

/**
 * Loads the Jones matrices onto the GPU.
 *
 * @param vm The VCSBeam context struct
 */
void vmPushJ( vcsbeam_context *vm )
{
/*#ifdef DEBUG
    printf("J_size_bytes = %lu\n", vm->J_size_bytes);
    uint8_t *dummy = (uint8_t *)(vm->J);
    uintptr_t byte;
    for (byte = 0; byte < vm->J_size_bytes; byte++)
    {
        if (byte % 16 == 0)  printf("\n");
        if (byte % 8 == 0)  printf(" ");
        printf ("%02x ", dummy[byte]);
    }
#endif */
    gpuMemcpy( vm->d_J, vm->J, vm->J_size_bytes, gpuMemcpyHostToDevice );
}

/**
 * Sets up CUDA streams for multi-pixel beamforming.
 *
 * @param vm The VCSBeam context struct
 *
 * \see vmDestroyCudaStreams()
 */
void vmCreateCudaStreams( vcsbeam_context *vm )
{
    vm->streams = (gpuStream_t *)malloc( vm->npointing * sizeof(gpuStream_t) );

    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        gpuStreamCreate( &(vm->streams[p]) );
    }
}

/**
 * Destroys the CUDA streams that were set up for multi-pixel beamforming.
 *
 * @param vm The VCSBeam context struct
 *
 * \see vmCreateCudaStreams()
 */
void vmDestroyCudaStreams( vcsbeam_context *vm )
{
    unsigned int p;
    for (p = 0; p < vm->npointing; p++)
    {
        gpuStreamDestroy( vm->streams[p] );
    }

    free( vm->streams );
}

/**
 * Allocates both CPU and GPU memory for the scales, offsets, and data for
 * PSRFITS output.
 *
 * @param vm The VCSBeam context struct
 * @param mpfs A pointer to a MPI PSRFITS struct
 *
 * \see vmDestroyStatistics()
 */
void vmCreateStatistics( vcsbeam_context *vm, mpi_psrfits *mpfs )
{
    uintptr_t nchan  = vm->nfine_chan;

    vm->offsets_size = vm->npointing*nchan*vm->out_nstokes*sizeof(float);
    vm->scales_size  = vm->npointing*nchan*vm->out_nstokes*sizeof(float);
    vm->Cscaled_size = vm->npointing*mpfs[0].coarse_chan_pf.sub.bytes_per_subint;

    gpuMalloc( (void **)&vm->d_offsets, vm->offsets_size );
    gpuMalloc( (void **)&vm->d_scales,  vm->scales_size );
    gpuMalloc( (void **)&vm->d_Cscaled, vm->Cscaled_size );

    gpuMallocHost( (void **)&vm->offsets, vm->offsets_size );
    gpuMallocHost( (void **)&vm->scales,  vm->scales_size );
    gpuMallocHost( (void **)&vm->Cscaled, vm->Cscaled_size );
}

/**
 * Frees both the CPU and GPU memory for the scales, offsets, and data for
 * PSRFITS output.
 *
 * @param vm The VCSBeam context struct
 *
 * \see vmCreateStatistics()
 */
void vmDestroyStatistics( vcsbeam_context *vm )
{
    gpuHostFree( vm->offsets );
    gpuHostFree( vm->scales );
    gpuHostFree( vm->Cscaled );

    gpuFree( vm->d_offsets );
    gpuFree( vm->d_scales );
    gpuFree( vm->d_Cscaled );
}

/**
 * Sets the number of pointings.
 *
 * @param vm The VCSBeam context struct
 * @param[in] npointings The number of pointings
 *
 * @todo Investigate whether this function is really needed
 */
void vmSetNumPointings( vcsbeam_context *vm, unsigned int npointings )
{
    uintptr_t npol      = vm->obs_metadata->num_ant_pols; // = 2
    vm->npointing       = npointings;
    vm->S_size_bytes    = vm->npointing * vm->nchan * vm->out_nstokes * vm->sample_rate * sizeof(float);
    vm->d_S_size_bytes  = vm->S_size_bytes;
    vm->e_size_bytes    = vm->npointing * vm->sample_rate * vm->nchan * npol * sizeof(gpuDoubleComplex);
    vm->d_e_size_bytes  = vm->e_size_bytes;
}

/**
 * Creates a list of file names for the input data.
 *
 * @param vm The VCSBeam context struct
 *
 * This function creates an array of file names that are passed to mwalib to
 * manage the reading in of the data.
 *
 * \see vmDestroyFilenames()
 */
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
        //printf("gps_second = file_start_timestep_second + t_idx*vm->seconds_per_file\n");

        gps_second = file_start_timestep_second + t_idx*vm->seconds_per_file;
        //printf("   %lu = %lu + %d*%d\n", gps_second, file_start_timestep_second, t_idx, vm->seconds_per_file);

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
 * Gets the input file name for the given channel index and GPS second.
 *
 * @param vm The VCSBeam context struct
 * @param[in] coarse_chan_idx The index of a coarse channel
 * @param[in] gps_second A GPS second
 * @param[out] filename A buffer for the filename
 *
 * This function returns the filename for the observation referred to in
 * `vm&rarr;obs_context`.
 *
 * The variable `vm&rarr;seconds_per_file` must be set to the relevant value
 * depending on whether the observation is Legacy (1) or MWAX (8), which is
 * done via vmLoadObsMetafits().
 *
 * `filename` must point to already-allocated memory.
 *
 * \see vmGetLegacyVoltFilename()
 */
void vmGetVoltFilename( vcsbeam_context *vm, unsigned int coarse_chan_idx, uint64_t gps_second, char *filename )
{
    uint64_t t0_gps_second = vm->obs_metadata->metafits_timesteps[0].gps_time_ms/1000;
    uintptr_t timestep_idx = (gps_second - t0_gps_second) / vm->seconds_per_file;

    // printf("vmGetVoltFilename: gps_second = %lu; t0_gps_second = %lu\n", gps_second, t0_gps_second);
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
 * Gets the file name for the given channel index and GPS, as if the
 * observation were a Legacy VCS observation.
 *
 * @param vm The VCSBeam context struct
 * @param[in] coarse_chan_idx The index of a coarse channel
 * @param[in] gps_second A GPS second
 * @param[out] filename A buffer for the filename
 *
 * This function returns the filename for the observation referred to in
 * `vm&rarr;obs_context_legacy`.
 *
 * `filename` must point to already-allocated memory.
 *
 * \see vmGetVoltFilename()
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

/**
 * Destroys the list of input file names.
 *
 * @param vm The VCSBeam context struct
 *
 * \see vmCreateFilenames()
 */
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

/**
 * Loads an observation's metadata from its metafits file.
 *
 * @param vm The VCSBeam context struct
 * @param filename The name of the metafits file to be loaded.
 *
 * This function loads the metadata into `vm&rarr;obs_context` and
 * `vm&rarr;obs_metadata`, using mwalib's API.
 * It also loads a "Legacy" version of the context and metafits into
 * `vm&rarr;obs_context_legacy` and `vm&rarr;obs_metafits_legacy`.
 * The observation to be loaded should be the "target" observation, i.e. the
 * observation whose VCS data is to be processed.
 */
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

    // Create OBS_METADATA_LEGACY
    if (mwalib_metafits_metadata_get( vm->obs_context_legacy, NULL, NULL, &vm->obs_metadata_legacy, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    vm->seconds_per_file = (vm->obs_metadata->mwa_version == VCSMWAXv2 ? 8 : 1);
                        // ^-- Raise an issue with mwalib (I shouldn't have to do this)

}

/**
 * Loads a calibration observation's metadata from its metafits file.
 *
 * @param vm The VCSBeam context struct
 * @param filename The name of the metafits file to be loaded.
 *
 * This function loads the metadata into `vm&rarr;cal_context` and
 * `vm&rarr;cal_metadata`, using mwalib's API.
 * The observation to be loaded should be the calibration observation, i.e.
 * the observation for which a calibration solution has been obtained.
 */
void vmLoadCalMetafits( vcsbeam_context *vm, char *filename )
{
    // Create CAL_CONTEXT
    if (mwalib_metafits_context_new2( filename, &vm->cal_context, vm->error_message, ERROR_MESSAGE_LEN) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits context: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }

    // Create CAL_METADATA
    if (mwalib_metafits_metadata_get( vm->cal_context, NULL, NULL, &vm->cal_metadata, vm->error_message, ERROR_MESSAGE_LEN ) != MWALIB_SUCCESS)
    {
        fprintf( stderr, "error (mwalib): cannot create metafits metadata: %s\n", vm->error_message );
        exit(EXIT_FAILURE);
    }
}

/**
 * Creates the voltage metadata structs using mwalib's API.
 *
 * @param vm The VCSBeam context struct
 *
 * This function should only be called after vmLoadObsMetafits().
 */
void vmGetVoltageMetadata( vcsbeam_context *vm )
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


/**
 * Gets the GPS second for an observation from a relative offset value.
 *
 * @param obs_metadata The observation's metadata
 * @param relative_begin An offset number of seconds
 * @return An absolute GPS second
 *
 * If `relative_begin` >= 0, then return the GPS second relative to the "good
 * time" (i.e. from the beginning of the observation).
 * If `relative_begin` < 0, then return the GPS second relative to the end of
 * the observation.
 *
 * | `relative_begin` | GPS second        |
 * | :--------------: | :---------------- |
 * | 0                | 1st "good" second |
 * | 1                | 2nd "good" second |
 * | 2                | 3rd "good" second |
 * | ...              | ...               |
 * | -2               | 2nd last second   |
 * | -1               | Last second       |
 */
long unsigned int get_relative_gps( MetafitsMetadata *obs_metadata, long int relative_begin )
{
    if (relative_begin >= 0)
        return obs_metadata->good_time_gps_ms/1000 + relative_begin;

    // Then relative begin must be negative, so count backwards from end
    return obs_metadata->metafits_timesteps[obs_metadata->num_metafits_timesteps + relative_begin].gps_time_ms/1000;
}

/**
 * Gets the GPS second from a string representation of either a relative or
 * absolute value.
 *
 * @param obs_metadata The observation's metadata
 * @param begin_str A string representation of either a relative or absolute
 *        GPS second
 * @return An absolute GPS second
 *
 * If the first character of `begin_str` is '`+`' or '`-`', then return a
 * relative GPS second according to get_relative_gps().
 * Otherwise, parse it as a GPS second in its own right.
 *
 * \see parse_coarse_chan_string()
 * \see get_relative_gps()
 */
long unsigned int parse_begin_string( MetafitsMetadata *obs_metadata, char *begin_str )
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


/**
 * Gets the coarse channel index from a string representation of either a relative or
 * absolute value.
 *
 * @param obs_metadata The observation's metadata
 * @param begin_coarse_chan_str A string representation of either a relative
 *        or absolute coarse channel index
 * @return An absolute coarse channel index
 *
 * If the first character of `begin_coarse_chan_str` is '`+`' or '`-`', then
 * return the coarse channel index relative to the lowest or highest coarse
 * channel, respectively (with "-1" representing the highest channel, "-2" the
 * second highest, etc.).
 * Otherwise, parse it as a coarse channel index in its own right.
 *
 * \see parse_begin_string()
 */
uintptr_t parse_coarse_chan_string( MetafitsMetadata *obs_metadata, char *begin_coarse_chan_str )
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

/**
 * Counts the number of tiles that are not flagged.
 *
 * @param vm The VCSBeam context struct
 *
 * The result is stored in `vm&rarr;num_not_flagged`.
 */
void vmSetNumNotFlaggedRFInputs( vcsbeam_context *vm )
{
    vm->num_not_flagged = 0;
    uintptr_t i;
    for (i = 0; i < vm->obs_metadata->num_rf_inputs; i++)
        if (!(vm->obs_metadata->rf_inputs[i].flagged))
            vm->num_not_flagged++;
}

/**
 * Finds a matching RF input in the given metadata.
 *
 * @param metadata The metadata to be searched
 * @param rfinput The RF input being sought
 * @return A pointer to the matching struct in `metadata`
 *
 * This function goes through the RF inputs in `metadata`, searching for one
 * that matches `rfinput`.
 * A "match" is an RF input that has the same `tile_id` and `pol`.
 * If no match is found, `NULL` is returned.
 */
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


/**
 * Finds a matching Antenna in the given metadata.
 *
 * @param metadata The metadata to be searched
 * @param rfinput The RF input being sought
 * @return A pointer to the matching struct in `metadata`
 *
 * This function goes through the antennas in `metadata`, searching for one
 * that matches `rfinput`.
 * A "match" is an RF input that has the same `tile_id`.
 * If no match is found, `NULL` is returned.
 */
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

/**
 * Gets the mwalib version.
 *
 * @param[out] version_str A buffer to hold the version string
 *
 * This function assumes that `version_str` is already allocated, and is big
 * enough
 */
void get_mwalib_version( char *version_str )
{
    sprintf( version_str, "%u.%u.%u",
            mwalib_get_version_major(),
            mwalib_get_version_minor(),
            mwalib_get_version_patch()
           );
}


/**
 * Parses RA/Dec pointings from a file.
 *
 * @param vm The VCSBeam context struct
 * @param filename The name of the file to be parsed.
 *
 * The file must contain whitespace-separated RAs and Decs in the format
 * `HH:MM:SS.S DD:MM:SS.S`.
 *
 * This function allocates memory for ras_hours and decs_degs arrays, which
 * will be destroyed during destroy_vcsbeam_context().
 */
void vmParsePointingFile( vcsbeam_context *vm, const char *filename )
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

/**
 * Reports all performance statistics.
 *
 * @param vm The VCSBeam context struct
 */
void vmReportPerformanceStats( vcsbeam_context *vm )
{
    logger_report_all_stats( vm->log );
}

/**
 * Prints a title to the specified log output stream.
 *
 * @param vm The VCSBeam context struct
 * @param title The text to be printed
 *
 * The text is printed in the format
 * ```
 * ------- VCSBeam (VERSION): TITLE -------
 * ```
 * with `VERSION` and `TITLE` being replaced with the VCSBeam version string
 * and the specified title respectively.
 */
void vmPrintTitle( vcsbeam_context *vm, const char *title )
{
    sprintf( vm->log_message, "------- VCSBeam (%s): %s -------",
            VCSBEAM_VERSION, title );
    logger_message( vm->log, vm->log_message );
}

void vmCheckError( vm_error err )
{
    switch (err)
    {
        case VM_END_OF_DATA:
            fprintf( stderr, "Error: End of data\n" );
            exit(EXIT_FAILURE);
            break;
        case VM_READ_BUFFER_NOT_SET:
            fprintf( stderr, "Error: Buffer not set\n" );
            exit(EXIT_FAILURE);
            break;
        case VM_READ_BUFFER_LOCKED:
            fprintf( stderr, "Error: Buffer locked\n" );
            exit(EXIT_FAILURE);
            break;
        default:
            break;
    }
}
