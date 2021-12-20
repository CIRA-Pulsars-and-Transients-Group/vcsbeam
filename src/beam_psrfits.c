/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <mwalib.h>
#include <star/pal.h>
#include <star/palmac.h>
#include <psrfits.h>

#include "vcsbeam.h"

void populate_spliced_psrfits_header(
        vcsbeam_context  *vm,
        struct psrfits   *pf,
        int               max_sec_per_file,
        int               outpol,
        beam_geom        *beam_geom_vals,
        char             *basename,
        bool              is_coherent )
{
    if ( !( outpol == 1 || outpol == 4 ) )
    {
        fprintf( stderr, "warning: populate_spliced_psrfits_header: "
                "unusual number of output pols = %d\n", outpol );
    }

    // Convert the UTC obs time into a string
    struct tm *ts = gmtime( &(vm->obs_metadata->sched_start_utc) );
    char   time_utc[64];
    strftime( time_utc, sizeof(time_utc), "%Y-%m-%dT%H:%M:%S", ts );

    // Get the sample rate
    unsigned int sample_rate = vm->fine_sample_rate;
    int coarse_chan_idx = vm->vcs_metadata->provided_coarse_chan_indices[0];
    int first_coarse_chan_idx = coarse_chan_idx - vm->mpi_rank;

    // Now set values for our hdrinfo structure
    strcpy( pf->hdr.project_id, vm->obs_metadata->project_id );
    strcpy( pf->hdr.obs_mode,  "SEARCH"    );
    strcpy( pf->hdr.observer,  "MWA User"  );
    strcpy( pf->hdr.telescope, "MWA"       );
    strcpy( pf->hdr.frontend,  "MWA-RECVR" );

    snprintf( pf->hdr.source, 24, "%u", vm->obs_metadata->obs_id );
    snprintf( pf->hdr.backend, 24, "%s", VCSBEAM_VERSION );

    pf->hdr.scanlen = 1.0; // in sec

    // Now let us finally get the time right
    strcpy(pf->hdr.date_obs,   time_utc);
    strcpy(pf->hdr.poln_type,  "LIN");
    strcpy(pf->hdr.track_mode, "TRACK");
    strcpy(pf->hdr.cal_mode,   "OFF");
    strcpy(pf->hdr.feed_mode,  "FA");

    pf->hdr.dt   = 1.0/sample_rate; // (sec)
    int last_coarse_chan_idx = first_coarse_chan_idx + vm->ncoarse_chans - 1;
    pf->hdr.fctr = 0.5*(
            vm->obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].chan_centre_hz +
            vm->obs_metadata->metafits_coarse_chans[last_coarse_chan_idx].chan_centre_hz) / 1e6; // (MHz)
    pf->hdr.BW = (
            vm->obs_metadata->metafits_coarse_chans[last_coarse_chan_idx].chan_end_hz -
            vm->obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].chan_start_hz
            ) / 1e6; // MHz

    // npols + nbits and whether pols are added
    pf->filenum       = 0;       // This is the crucial one to set to initialize things
    pf->rows_per_file = max_sec_per_file;     // I assume this is a max subint issue

    pf->hdr.npol         = outpol;
    pf->hdr.nchan        = vm->nfine_chan * vm->ncoarse_chans;
    pf->hdr.onlyI        = 0;

    pf->hdr.scan_number   = 1;
    pf->hdr.rcvr_polns    = 2;
    pf->hdr.offset_subint = 0;

    if (is_coherent)
        pf->hdr.summed_polns = 0;
    else
        pf->hdr.summed_polns = 1;

    uint32_t fine_chan_width = vm->obs_metadata->coarse_chan_width_hz / vm->nfine_chan;
    pf->hdr.df         = fine_chan_width / 1e6; // (MHz)
    pf->hdr.orig_nchan = pf->hdr.nchan;
    pf->hdr.orig_df    = pf->hdr.df;
    pf->hdr.nbits      = 8;
    pf->hdr.nsblk      = sample_rate;  // block is always 1 second of data

    pf->hdr.ds_freq_fact = 1;
    pf->hdr.ds_time_fact = 1;

    // some things that we are unlikely to change
    pf->hdr.fd_hand  = 1;
    pf->hdr.fd_sang  = 45.0;
    pf->hdr.fd_xyph  = 0.0;
    pf->hdr.be_phase = 0;
    pf->hdr.chan_dm  = 0.0;

    // Now set values for our subint structure
    pf->tot_rows     = 0;
    pf->sub.tsubint  = roundf(pf->hdr.nsblk * pf->hdr.dt);
    pf->sub.offs     = roundf(pf->tot_rows * pf->sub.tsubint) + 0.5*pf->sub.tsubint;

    pf->sub.feed_ang = 0.0;
    pf->sub.pos_ang  = 0.0;
    pf->sub.par_ang  = 0.0;

    // Specify psrfits data type
    pf->sub.FITS_typecode = TBYTE;

    pf->sub.bytes_per_subint = (pf->hdr.nbits * pf->hdr.nchan *
            pf->hdr.npol  * pf->hdr.nsblk) / 8;

    // Create and initialize the subint arrays
    pf->sub.dat_freqs   = (float *)malloc(sizeof(float) * pf->hdr.nchan);
    pf->sub.dat_weights = (float *)malloc(sizeof(float) * pf->hdr.nchan);

    int i; // idx into dat_freqs
    int iF, iC; // mwalib (i)dxs for (F)ine and (C)oarse channels
    uint32_t start_hz = vm->obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].chan_start_hz;
    for (i = 0 ; i < pf->hdr.nchan; i++)
    {
        iC = i / vm->nfine_chan + first_coarse_chan_idx;
        iF = (iC * vm->nfine_chan) + (i % vm->nfine_chan);
        pf->sub.dat_freqs[i] = start_hz + iF*fine_chan_width;
        pf->sub.dat_weights[i] = 1.0;
    }

    // the following is definitely needed for 8 bit numbers
    pf->sub.dat_offsets = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);
    pf->sub.dat_scales  = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);

    pf->sub.data    = (unsigned char *)malloc(pf->sub.bytes_per_subint);
    pf->sub.rawdata = pf->sub.data;

    // Update values that depend on get_jones()
    if (beam_geom_vals != NULL)
    {
        if (is_coherent) 
        {
            pf->hdr.ra2000  = beam_geom_vals->mean_ra  * R2D;
            pf->hdr.dec2000 = beam_geom_vals->mean_dec * R2D;
        } 
        else
        {
            // Use the tile pointing instead of the pencil beam pointing
            pf->hdr.ra2000  = vm->obs_metadata->ra_tile_pointing_deg;
            pf->hdr.dec2000 = vm->obs_metadata->dec_tile_pointing_deg;
        }

        dec2hms( pf->hdr.ra_str,  pf->hdr.ra2000/15.0, 0 );
        dec2hms( pf->hdr.dec_str, pf->hdr.dec2000,     1 );

        pf->hdr.azimuth    = beam_geom_vals->az*R2D;
        pf->hdr.zenith_ang = 90.0 - (beam_geom_vals->el*R2D);

        pf->hdr.beam_FWHM = 0.25;
        pf->hdr.start_lst = beam_geom_vals->lmst * 60.0 * 60.0;        // Local Apparent Sidereal Time in seconds
        pf->hdr.start_sec = roundf(beam_geom_vals->fracmjd*86400.0);   // this will always be a whole second
        pf->hdr.start_day = beam_geom_vals->intmjd;
        pf->hdr.MJD_epoch = beam_geom_vals->intmjd + beam_geom_vals->fracmjd;

        // Now set values for our subint structure
        pf->sub.lst      = pf->hdr.start_lst;
        pf->sub.ra       = pf->hdr.ra2000;
        pf->sub.dec      = pf->hdr.dec2000;
        palEqgal(pf->hdr.ra2000*D2R, pf->hdr.dec2000*D2R,
                &pf->sub.glon, &pf->sub.glat);
        pf->sub.glon    *= R2D;
        pf->sub.glat    *= R2D;
        pf->sub.tel_az   = pf->hdr.azimuth;
        pf->sub.tel_zen  = pf->hdr.zenith_ang;

        // Construct the output filename
        if (basename != NULL)
            strcpy( pf->basefilename, basename );
        else
        {
            char chan_str[16];
            if (first_coarse_chan_idx == last_coarse_chan_idx)
            {
                sprintf( chan_str, "%03ld",
                        vm->obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].rec_chan_number );
            }
            else
            {
                sprintf( chan_str, "%03ld-%03ld",
                        vm->obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].rec_chan_number,
                        vm->obs_metadata->metafits_coarse_chans[last_coarse_chan_idx].rec_chan_number );
            }

            if (is_coherent)
            {
                sprintf( pf->basefilename, "%s_%s_%s_%s_ch%s",
                        pf->hdr.project_id,
                        pf->hdr.source, 
                        pf->hdr.ra_str, pf->hdr.dec_str,
                        chan_str );
            }
            else
            {
                sprintf( pf->basefilename, "%s_%s_ch%s_incoh",
                        pf->hdr.project_id,
                        pf->hdr.source,
                        chan_str );
            }
        }
    }
}

void populate_psrfits_header(
        vcsbeam_context  *vm,
        struct psrfits   *pf,
        int               max_sec_per_file,
        int               outpol,
        beam_geom        *beam_geom_vals,
        char             *incoh_basename,
        bool              is_coherent )
{
    if ( !( outpol == 1 || outpol == 4 ) )
    {
        fprintf( stderr, "warning: populate_psrfits_header: "
                "unusual number of output pols = %d\n", outpol );
    }

    // Convert the UTC obs time into a string
    struct tm *ts = gmtime( &(vm->obs_metadata->sched_start_utc) );
    char   time_utc[64];
    strftime( time_utc, sizeof(time_utc), "%Y-%m-%dT%H:%M:%S", ts );

    // Get the sample rate
    unsigned int sample_rate = vm->fine_sample_rate;
    int coarse_chan_idx = vm->vcs_metadata->provided_coarse_chan_indices[0];

    // Now set values for our hdrinfo structure
    strcpy( pf->hdr.project_id, vm->obs_metadata->project_id );
    strcpy( pf->hdr.obs_mode,  "SEARCH"    );
    strcpy( pf->hdr.observer,  "MWA User"  );
    strcpy( pf->hdr.telescope, "MWA"       );
    strcpy( pf->hdr.frontend,  "MWA-RECVR" );

    snprintf( pf->hdr.source, 24, "%u", vm->obs_metadata->obs_id );
    snprintf( pf->hdr.backend, 24, "%s", VCSBEAM_VERSION );

    pf->hdr.scanlen = 1.0; // in sec

    // Now let us finally get the time right
    strcpy(pf->hdr.date_obs,   time_utc);
    strcpy(pf->hdr.poln_type,  "LIN");
    strcpy(pf->hdr.track_mode, "TRACK");
    strcpy(pf->hdr.cal_mode,   "OFF");
    strcpy(pf->hdr.feed_mode,  "FA");

    pf->hdr.dt   = 1.0/sample_rate; // (sec)
    pf->hdr.fctr = vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_centre_hz / 1e6; // (MHz)
    pf->hdr.BW   = vm->obs_metadata->coarse_chan_width_hz / 1e6;  // (MHz)

    // npols + nbits and whether pols are added
    pf->filenum       = 0;       // This is the crucial one to set to initialize things
    pf->rows_per_file = max_sec_per_file;     // I assume this is a max subint issue

    pf->hdr.npol         = outpol;
    pf->hdr.nchan        = vm->nfine_chan;
    pf->hdr.onlyI        = 0;

    pf->hdr.scan_number   = 1;
    pf->hdr.rcvr_polns    = 2;
    pf->hdr.offset_subint = 0;

    if (is_coherent)
        pf->hdr.summed_polns = 0;
    else
        pf->hdr.summed_polns = 1;

    uint32_t fine_chan_width = vm->obs_metadata->coarse_chan_width_hz / vm->nfine_chan;
    pf->hdr.df         = fine_chan_width / 1e6; // (MHz)
    pf->hdr.orig_nchan = pf->hdr.nchan;
    pf->hdr.orig_df    = pf->hdr.df;
    pf->hdr.nbits      = 8;
    pf->hdr.nsblk      = sample_rate;  // block is always 1 second of data

    pf->hdr.ds_freq_fact = 1;
    pf->hdr.ds_time_fact = 1;

    // some things that we are unlikely to change
    pf->hdr.fd_hand  = 1;
    pf->hdr.fd_sang  = 45.0;
    pf->hdr.fd_xyph  = 0.0;
    pf->hdr.be_phase = 0;
    pf->hdr.chan_dm  = 0.0;

    // Now set values for our subint structure
    pf->tot_rows     = 0;
    pf->sub.tsubint  = roundf(pf->hdr.nsblk * pf->hdr.dt);
    pf->sub.offs     = roundf(pf->tot_rows * pf->sub.tsubint) + 0.5*pf->sub.tsubint;

    pf->sub.feed_ang = 0.0;
    pf->sub.pos_ang  = 0.0;
    pf->sub.par_ang  = 0.0;

    // Specify psrfits data type
    pf->sub.FITS_typecode = TBYTE;

    pf->sub.bytes_per_subint = (pf->hdr.nbits * pf->hdr.nchan *
            pf->hdr.npol  * pf->hdr.nsblk) / 8;

    // Create and initialize the subint arrays
    pf->sub.dat_freqs   = (float *)malloc(sizeof(float) * pf->hdr.nchan);
    pf->sub.dat_weights = (float *)malloc(sizeof(float) * pf->hdr.nchan);

    double dtmp = pf->hdr.fctr - 0.5 * pf->hdr.BW + 0.5 * pf->hdr.df;
    int i;
    for (i = 0 ; i < pf->hdr.nchan ; i++) {
        pf->sub.dat_freqs[i] = dtmp + i * pf->hdr.df;
        pf->sub.dat_weights[i] = 1.0;
    }

    // the following is definitely needed for 8 bit numbers
    pf->sub.dat_offsets = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);
    pf->sub.dat_scales  = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);
    for (i = 0 ; i < pf->hdr.nchan * pf->hdr.npol ; i++) {
        pf->sub.dat_offsets[i] = 0.0;
        pf->sub.dat_scales[i]  = 1.0;
    }

    pf->sub.data    = (unsigned char *)malloc(pf->sub.bytes_per_subint);
    pf->sub.rawdata = pf->sub.data;

    // Update values that depend on get_jones()
    if (beam_geom_vals != NULL)
    {
        if (is_coherent) 
        {
            pf->hdr.ra2000  = beam_geom_vals->mean_ra  * R2D;
            pf->hdr.dec2000 = beam_geom_vals->mean_dec * R2D;
        } 
        else
        {
            // Use the tile pointing instead of the pencil beam pointing
            pf->hdr.ra2000  = vm->obs_metadata->ra_tile_pointing_deg;
            pf->hdr.dec2000 = vm->obs_metadata->dec_tile_pointing_deg;
        }

        dec2hms( pf->hdr.ra_str,  pf->hdr.ra2000/15.0, 0 );
        dec2hms( pf->hdr.dec_str, pf->hdr.dec2000,     1 );

        pf->hdr.azimuth    = beam_geom_vals->az*R2D;
        pf->hdr.zenith_ang = 90.0 - (beam_geom_vals->el*R2D);

        pf->hdr.beam_FWHM = 0.25;
        pf->hdr.start_lst = beam_geom_vals->lmst * 60.0 * 60.0;        // Local Apparent Sidereal Time in seconds
        pf->hdr.start_sec = roundf(beam_geom_vals->fracmjd*86400.0);   // this will always be a whole second
        pf->hdr.start_day = beam_geom_vals->intmjd;
        pf->hdr.MJD_epoch = beam_geom_vals->intmjd + beam_geom_vals->fracmjd;

        // Now set values for our subint structure
        pf->sub.lst      = pf->hdr.start_lst;
        pf->sub.ra       = pf->hdr.ra2000;
        pf->sub.dec      = pf->hdr.dec2000;
        palEqgal(pf->hdr.ra2000*D2R, pf->hdr.dec2000*D2R,
                &pf->sub.glon, &pf->sub.glat);
        pf->sub.glon    *= R2D;
        pf->sub.glat    *= R2D;
        pf->sub.tel_az   = pf->hdr.azimuth;
        pf->sub.tel_zen  = pf->hdr.zenith_ang;

        if (is_coherent)
        {
            sprintf( pf->basefilename, "%s_%s_%s_%s_ch%03ld",
                    pf->hdr.project_id,
                    pf->hdr.source, 
                    pf->hdr.ra_str, pf->hdr.dec_str,
                    vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
        }
        else
        {
            if (incoh_basename != NULL)
                strcpy( pf->basefilename, incoh_basename );
            else
                sprintf( pf->basefilename, "%s_%s_incoh_ch%03ld",
                        pf->hdr.project_id,
                        pf->hdr.source,
                        vm->obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
        }
    }
}


void free_psrfits( struct psrfits *pf )
{
    free( pf->sub.data        );
    free( pf->sub.dat_freqs   );
    free( pf->sub.dat_weights );
    free( pf->sub.dat_offsets );
    free( pf->sub.dat_scales  );
}


float *create_data_buffer_psrfits( size_t size )
{
    float *ptr = (float *)malloc( size * sizeof(float) );
    return ptr;
}

void vmInitMPIPsrfits(
        vcsbeam_context *vm,
        mpi_psrfits *mpf,
        int max_sec_per_file,
        int nstokes,
        beam_geom *bg,
        char *outfile,
        bool is_coherent )
{
    if (vm->mpi_rank == PERFORMANCE_NO_MPI)
    {
        fprintf( stderr, "error: init_mpi_psrfits: MPI not being used\n" );
        exit(EXIT_FAILURE);
    }

    mpf->writer_id = vm->writer;

    // Populate the PSRFITS header struct for the combined (spliced) output file
    if (vm->mpi_rank == mpf->writer_id)
        populate_spliced_psrfits_header( vm, &(mpf->spliced_pf),
                max_sec_per_file, nstokes,
                bg, outfile, is_coherent );

    // Populate the PSRFITS header struct for a single channel
    populate_psrfits_header( vm, &(mpf->coarse_chan_pf), max_sec_per_file,
            nstokes, bg, outfile, is_coherent );

    // Create MPI vector types designed to splice the coarse channels together
    // correctly during MPI_Gather
    MPI_Type_contiguous( mpf->coarse_chan_pf.hdr.nchan, MPI_BYTE, &(mpf->coarse_chan_spectrum) );
    MPI_Type_commit( &(mpf->coarse_chan_spectrum) );

    MPI_Type_vector( mpf->coarse_chan_pf.hdr.nsblk*nstokes, 1, vm->ncoarse_chans,
            mpf->coarse_chan_spectrum, &(mpf->total_spectrum_type) );
    MPI_Type_commit( &(mpf->total_spectrum_type) );

    MPI_Type_create_resized( mpf->total_spectrum_type, 0, mpf->coarse_chan_pf.hdr.nchan, &(mpf->spliced_type) );
    MPI_Type_commit( &(mpf->spliced_type) );
}

void free_mpi_psrfits( mpi_psrfits *mpf )
{
    MPI_Type_free( &(mpf->spliced_type) );
    MPI_Type_free( &(mpf->total_spectrum_type) );
    MPI_Type_free( &(mpf->coarse_chan_spectrum) );

    free_psrfits( &(mpf->coarse_chan_pf) );

    int mpi_proc_id;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_proc_id );
    if (mpi_proc_id == mpf->writer_id)
    {
        free_psrfits( &(mpf->spliced_pf) );
    }
}

void gather_splice_psrfits( mpi_psrfits *mpf )
{
    int nsamples = mpf->coarse_chan_pf.hdr.nsblk;
    int nstokes = mpf->coarse_chan_pf.hdr.npol;
    int nfinechans = mpf->coarse_chan_pf.hdr.nchan;

    MPI_Igather(
            mpf->coarse_chan_pf.sub.data, nsamples*nstokes, mpf->coarse_chan_spectrum,
            mpf->spliced_pf.sub.data, 1, mpf->spliced_type,
            mpf->writer_id, MPI_COMM_WORLD, &(mpf->request_data) );

    MPI_Igather(
            mpf->coarse_chan_pf.sub.dat_offsets, nfinechans*nstokes, MPI_FLOAT,
            mpf->spliced_pf.sub.dat_offsets, nfinechans*nstokes, MPI_FLOAT,
            mpf->writer_id, MPI_COMM_WORLD, &(mpf->request_offsets) );

    MPI_Igather(
            mpf->coarse_chan_pf.sub.dat_scales, nfinechans*nstokes, MPI_FLOAT,
            mpf->spliced_pf.sub.dat_scales, nfinechans*nstokes, MPI_FLOAT,
            mpf->writer_id, MPI_COMM_WORLD, &(mpf->request_scales) );
}

void wait_splice_psrfits( mpi_psrfits *mpf )
{
    MPI_Wait( &(mpf->request_data),    MPI_STATUS_IGNORE );
    MPI_Wait( &(mpf->request_offsets), MPI_STATUS_IGNORE );
    MPI_Wait( &(mpf->request_scales),  MPI_STATUS_IGNORE );
}
