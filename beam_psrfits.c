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
#include <mwalib.h>
#include <mpi.h>
#include "star/pal.h"
#include "star/palmac.h"
#include "psrfits.h"
#include "beam_psrfits.h"
#include "geometry.h"
#include "metadata.h"
#include "jones.h"

void float_to_unit8(float * in, int n, int8_t *out)
{
    int j;
    float min = -128.0; // -126.0 and -128.0 give the same result on test data
    float max = 127.0;
    // use a temp var so we don't modify the input data
    float scratch;
    for (j = 0; j < n; j++) {
        // TODO: count the number of samples that were clipped, store that and put it in the psrfits header
        // the if branching and ternary updates seem to be equivalent execution time
        if (in[j]> max) {
            scratch = max;
        } else if (in[j] < min) {
            scratch = min;
        } else {
            scratch = in[j];
        }
//        scratch = (in[j] > max) ? (max) : in[j];
//        scratch = (in[j] < min) ? (min) : scratch;
        out[j] = (uint8_t)( (int8_t)rint(scratch) + 128);
    }

}


void populate_spliced_psrfits_header(
        struct psrfits   *pf,
        MetafitsMetadata *obs_metadata,
        VoltageMetadata  *vcs_metadata,
        int               first_coarse_chan_idx,
        int               ncoarse_chans,
        int               max_sec_per_file,
        int               outpol,
        struct beam_geom *beam_geom_vals,
        char             *basename,
        bool              is_coherent )
{
    if ( !( outpol == 1 || outpol == 4 ) )
    {
        fprintf( stderr, "warning: populate_psrfits_header: "
                "unusual number of output pols = %d\n", outpol );
    }

    // Convert the UTC obs time into a string
    struct tm *ts = gmtime( &(obs_metadata->sched_start_utc) );
    char   time_utc[64];
    strftime( time_utc, sizeof(time_utc), "%Y-%m-%dT%H:%M:%S", ts );

    // Get the sample rate
    unsigned int sample_rate = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;

    // Now set values for our hdrinfo structure
    strcpy( pf->hdr.project_id, obs_metadata->project_id );
    strcpy( pf->hdr.obs_mode,  "SEARCH"    );
    strcpy( pf->hdr.observer,  "MWA User"  );
    strcpy( pf->hdr.telescope, "MWA"       );
    strcpy( pf->hdr.frontend,  "MWA-RECVR" );

    snprintf( pf->hdr.source, 24, "%u", obs_metadata->obs_id );
    snprintf( pf->hdr.backend, 24, "%s", VERSION_BEAMFORMER );

    pf->hdr.scanlen = 1.0; // in sec

    // Now let us finally get the time right
    strcpy(pf->hdr.date_obs,   time_utc);
    strcpy(pf->hdr.poln_type,  "LIN");
    strcpy(pf->hdr.track_mode, "TRACK");
    strcpy(pf->hdr.cal_mode,   "OFF");
    strcpy(pf->hdr.feed_mode,  "FA");

    pf->hdr.dt   = 1.0/sample_rate; // (sec)
    int last_coarse_chan_idx = first_coarse_chan_idx + ncoarse_chans - 1;
    pf->hdr.fctr = 0.5*(
            obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].chan_centre_hz +
            obs_metadata->metafits_coarse_chans[last_coarse_chan_idx].chan_centre_hz) / 1e6; // (MHz)
    pf->hdr.BW = (
            obs_metadata->metafits_coarse_chans[last_coarse_chan_idx].chan_end_hz -
            obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].chan_start_hz
            ) / 1e6; // MHz

    // npols + nbits and whether pols are added
    pf->filenum       = 0;       // This is the crucial one to set to initialize things
    pf->rows_per_file = max_sec_per_file;     // I assume this is a max subint issue

    pf->hdr.npol         = outpol;
    pf->hdr.nchan        = obs_metadata->num_volt_fine_chans_per_coarse * ncoarse_chans;
    pf->hdr.onlyI        = 0;

    pf->hdr.scan_number   = 1;
    pf->hdr.rcvr_polns    = 2;
    pf->hdr.offset_subint = 0;

    if (is_coherent)
        pf->hdr.summed_polns = 0;
    else
        pf->hdr.summed_polns = 1;

    pf->hdr.df         = obs_metadata->volt_fine_chan_width_hz / 1e6; // (MHz)
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
    for (i = 0 ; i < pf->hdr.nchan; i++)
    {
        iC = i / obs_metadata->num_volt_fine_chans_per_coarse + first_coarse_chan_idx;
        iF = iC*obs_metadata->num_volt_fine_chans_per_coarse +
            i % obs_metadata->num_volt_fine_chans_per_coarse;
        pf->sub.dat_freqs[i] = obs_metadata->metafits_fine_chan_freqs_hz[iF] / 1e6;
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
            pf->hdr.ra2000  = beam_geom_vals->mean_ra  * PAL__DR2D;
            pf->hdr.dec2000 = beam_geom_vals->mean_dec * PAL__DR2D;
        } 
        else
        {
            // Use the tile pointing instead of the pencil beam pointing
            pf->hdr.ra2000  = obs_metadata->ra_tile_pointing_deg;
            pf->hdr.dec2000 = obs_metadata->dec_tile_pointing_deg;
        }

        dec2hms( pf->hdr.ra_str,  pf->hdr.ra2000/15.0, 0 );
        dec2hms( pf->hdr.dec_str, pf->hdr.dec2000,     1 );

        pf->hdr.azimuth    = beam_geom_vals->az*PAL__DR2D;
        pf->hdr.zenith_ang = 90.0 - (beam_geom_vals->el*PAL__DR2D);

        pf->hdr.beam_FWHM = 0.25;
        pf->hdr.start_lst = beam_geom_vals->lmst * 60.0 * 60.0;        // Local Apparent Sidereal Time in seconds
        pf->hdr.start_sec = roundf(beam_geom_vals->fracmjd*86400.0);   // this will always be a whole second
        pf->hdr.start_day = beam_geom_vals->intmjd;
        pf->hdr.MJD_epoch = beam_geom_vals->intmjd + beam_geom_vals->fracmjd;

        // Now set values for our subint structure
        pf->sub.lst      = pf->hdr.start_lst;
        pf->sub.ra       = pf->hdr.ra2000;
        pf->sub.dec      = pf->hdr.dec2000;
        palEqgal(pf->hdr.ra2000*PAL__DD2R, pf->hdr.dec2000*PAL__DD2R,
                &pf->sub.glon, &pf->sub.glat);
        pf->sub.glon    *= PAL__DR2D;
        pf->sub.glat    *= PAL__DR2D;
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
                        obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].rec_chan_number );
            }
            else
            {
                sprintf( chan_str, "%03ld-%03ld",
                        obs_metadata->metafits_coarse_chans[first_coarse_chan_idx].rec_chan_number,
                        obs_metadata->metafits_coarse_chans[last_coarse_chan_idx].rec_chan_number );
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
        struct psrfits   *pf,
        MetafitsMetadata *obs_metadata,
        VoltageMetadata  *vcs_metadata,
        int               coarse_chan_idx,
        int               max_sec_per_file,
        int               outpol,
        struct beam_geom *beam_geom_vals,
        char             *incoh_basename,
        bool              is_coherent )
{
    if ( !( outpol == 1 || outpol == 4 ) )
    {
        fprintf( stderr, "warning: populate_psrfits_header: "
                "unusual number of output pols = %d\n", outpol );
    }

    // Convert the UTC obs time into a string
    struct tm *ts = gmtime( &(obs_metadata->sched_start_utc) );
    char   time_utc[64];
    strftime( time_utc, sizeof(time_utc), "%Y-%m-%dT%H:%M:%S", ts );

    // Get the sample rate
    unsigned int sample_rate = vcs_metadata->num_samples_per_voltage_block * vcs_metadata->num_voltage_blocks_per_second;

    // Now set values for our hdrinfo structure
    strcpy( pf->hdr.project_id, obs_metadata->project_id );
    strcpy( pf->hdr.obs_mode,  "SEARCH"    );
    strcpy( pf->hdr.observer,  "MWA User"  );
    strcpy( pf->hdr.telescope, "MWA"       );
    strcpy( pf->hdr.frontend,  "MWA-RECVR" );

    snprintf( pf->hdr.source, 24, "%u", obs_metadata->obs_id );
    snprintf( pf->hdr.backend, 24, "%s", VERSION_BEAMFORMER );

    pf->hdr.scanlen = 1.0; // in sec

    // Now let us finally get the time right
    strcpy(pf->hdr.date_obs,   time_utc);
    strcpy(pf->hdr.poln_type,  "LIN");
    strcpy(pf->hdr.track_mode, "TRACK");
    strcpy(pf->hdr.cal_mode,   "OFF");
    strcpy(pf->hdr.feed_mode,  "FA");

    pf->hdr.dt   = 1.0/sample_rate; // (sec)
    pf->hdr.fctr = obs_metadata->metafits_coarse_chans[coarse_chan_idx].chan_centre_hz / 1e6; // (MHz)
    pf->hdr.BW   = obs_metadata->coarse_chan_width_hz / 1e6;  // (MHz)

    // npols + nbits and whether pols are added
    pf->filenum       = 0;       // This is the crucial one to set to initialize things
    pf->rows_per_file = max_sec_per_file;     // I assume this is a max subint issue

    pf->hdr.npol         = outpol;
    pf->hdr.nchan        = obs_metadata->num_volt_fine_chans_per_coarse;
    pf->hdr.onlyI        = 0;

    pf->hdr.scan_number   = 1;
    pf->hdr.rcvr_polns    = 2;
    pf->hdr.offset_subint = 0;

    if (is_coherent)
        pf->hdr.summed_polns = 0;
    else
        pf->hdr.summed_polns = 1;

    pf->hdr.df         = obs_metadata->volt_fine_chan_width_hz / 1e6; // (MHz)
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
            pf->hdr.ra2000  = beam_geom_vals->mean_ra  * PAL__DR2D;
            pf->hdr.dec2000 = beam_geom_vals->mean_dec * PAL__DR2D;
        } 
        else
        {
            // Use the tile pointing instead of the pencil beam pointing
            pf->hdr.ra2000  = obs_metadata->ra_tile_pointing_deg;
            pf->hdr.dec2000 = obs_metadata->dec_tile_pointing_deg;
        }

        dec2hms( pf->hdr.ra_str,  pf->hdr.ra2000/15.0, 0 );
        dec2hms( pf->hdr.dec_str, pf->hdr.dec2000,     1 );

        pf->hdr.azimuth    = beam_geom_vals->az*PAL__DR2D;
        pf->hdr.zenith_ang = 90.0 - (beam_geom_vals->el*PAL__DR2D);

        pf->hdr.beam_FWHM = 0.25;
        pf->hdr.start_lst = beam_geom_vals->lmst * 60.0 * 60.0;        // Local Apparent Sidereal Time in seconds
        pf->hdr.start_sec = roundf(beam_geom_vals->fracmjd*86400.0);   // this will always be a whole second
        pf->hdr.start_day = beam_geom_vals->intmjd;
        pf->hdr.MJD_epoch = beam_geom_vals->intmjd + beam_geom_vals->fracmjd;

        // Now set values for our subint structure
        pf->sub.lst      = pf->hdr.start_lst;
        pf->sub.ra       = pf->hdr.ra2000;
        pf->sub.dec      = pf->hdr.dec2000;
        palEqgal(pf->hdr.ra2000*PAL__DD2R, pf->hdr.dec2000*PAL__DD2R,
                &pf->sub.glon, &pf->sub.glat);
        pf->sub.glon    *= PAL__DR2D;
        pf->sub.glat    *= PAL__DR2D;
        pf->sub.tel_az   = pf->hdr.azimuth;
        pf->sub.tel_zen  = pf->hdr.zenith_ang;

        if (is_coherent)
        {
            sprintf( pf->basefilename, "%s_%s_%s_%s_ch%03ld",
                    pf->hdr.project_id,
                    pf->hdr.source, 
                    pf->hdr.ra_str, pf->hdr.dec_str,
                    obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
        }
        else
        {
            if (incoh_basename != NULL)
                strcpy( pf->basefilename, incoh_basename );
            else
                sprintf( pf->basefilename, "%s_%s_incoh_ch%03ld",
                        pf->hdr.project_id,
                        pf->hdr.source,
                        obs_metadata->metafits_coarse_chans[coarse_chan_idx].rec_chan_number );
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

void correct_psrfits_stt( struct psrfits *pf )
{
    /* now we have to correct the STT_SMJD/STT_OFFS as they will have been broken by the write_psrfits*/
    int    itmp    = 0;
    int    itmp2   = 0;
    double dtmp    = 0;
    int    status  = 0;

    //fits_open_file(&(pf.fptr), pf.filename, READWRITE, &status);

    fits_read_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp,  NULL, &status);
    fits_read_key(pf->fptr, TINT,    "STT_SMJD", &itmp,  NULL, &status);
    fits_read_key(pf->fptr, TINT,    "STT_IMJD", &itmp2, NULL, &status);

    if (dtmp > 0.5) {
        itmp = itmp+1;
        if (itmp == 86400) {
            itmp = 0;
            itmp2++;
        }
    }
    dtmp = 0.0;

    fits_update_key(pf->fptr, TINT, "STT_SMJD", &itmp, NULL, &status);
    fits_update_key(pf->fptr, TINT, "STT_IMJD", &itmp2, NULL, &status);
    fits_update_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp, NULL, &status);

}


void psrfits_write_second( struct psrfits *pf, float *data_buffer, int nchan,
                           int outpol, int p )
{
    int sec_size = outpol * nchan * pf->hdr.nsblk;
    int8_t *out_buffer_8 = (int8_t *)malloc( sec_size * sizeof(int8_t) );

    // pointing_offset makes the buffer start at the memory assigned the pointing
    int pointing_offset = p * sec_size;
    float *pointing_buffer  = malloc( sec_size * sizeof(float) );
    memcpy(pointing_buffer, data_buffer + pointing_offset, sec_size * sizeof(float) );
    float_to_unit8( pointing_buffer, sec_size, out_buffer_8);
    
    memcpy( pf->sub.data, out_buffer_8, pf->sub.bytes_per_subint );
    //memset(pf->filename,0,strlen(pf->filename));
    //memset(pf->hdr.poln_order,0,strlen(pf->hdr.poln_order));
    
    if (psrfits_write_subint(pf) != 0)
    {
        fprintf(stderr, "error: Write subint failed. File exists?\n");
        exit(EXIT_FAILURE);
    }
    pf->sub.offs = roundf(pf->tot_rows * pf->sub.tsubint) + 0.5*pf->sub.tsubint;
    pf->sub.lst += pf->sub.tsubint;
    
    free( pointing_buffer );
    free( out_buffer_8 );
}


void printf_psrfits( struct psrfits *pf )
{
    fprintf(stdout, "\nPSRFITS:\n");
    fprintf(stdout, "Basename of output file     [%s]\n", pf->basefilename);
    fprintf(stdout, "Filename of output file     [%s]\n", pf->filename);
    fprintf(stdout, "CFITSIO file pointer        [%p]\n", pf->fptr);

    fprintf(stdout, "\nPSRFITS HDRINFO:\n");
    fprintf(stdout, "Obs mode                    [%s]\n", pf->hdr.obs_mode);
    fprintf(stdout, "Telescope                   [%s]\n", pf->hdr.telescope);
    fprintf(stdout, "Observer                    [%s]\n", pf->hdr.observer);
    fprintf(stdout, "Source                      [%s]\n", pf->hdr.source);
    fprintf(stdout, "Front End                   [%s]\n", pf->hdr.frontend);
    fprintf(stdout, "Back End                    [%s]\n", pf->hdr.backend);
    fprintf(stdout, "Project ID                  [%s]\n", pf->hdr.project_id);
    fprintf(stdout, "Date Obs                    [%s]\n", pf->hdr.date_obs);
    fprintf(stdout, "RA (string)                 [%s]\n", pf->hdr.ra_str);
    fprintf(stdout, "Dec (string)                [%s]\n", pf->hdr.dec_str);
    fprintf(stdout, "Pol recorded (LIN or CIRC)  [%s]\n", pf->hdr.poln_type);
    fprintf(stdout, "Order of pols               [%s]\n", pf->hdr.poln_order);
    fprintf(stdout, "Track mode                  [%s]\n", pf->hdr.track_mode);
    fprintf(stdout, "Cal mode                    [%s]\n", pf->hdr.cal_mode);
    fprintf(stdout, "Feed mode                   [%s]\n", pf->hdr.feed_mode);
    fprintf(stdout, "Start MJD                   [%Lf]\n", pf->hdr.MJD_epoch);
    fprintf(stdout, "Sample Time (s)             [%lf]\n", pf->hdr.dt);
    fprintf(stdout, "Centre Frequency (MHz)      [%lf]\n", pf->hdr.fctr);
    fprintf(stdout, "Orig freq spacing (MHz)     [%lf]\n", pf->hdr.orig_df);
    fprintf(stdout, "Freq spacing (MHz)          [%lf]\n", pf->hdr.df);
    fprintf(stdout, "Bandwidth (MHz)             [%lf]\n", pf->hdr.BW);
    fprintf(stdout, "RA (2000) (deg)             [%lf]\n", pf->hdr.ra2000);
    fprintf(stdout, "Dec (2000) (deg)            [%lf]\n", pf->hdr.dec2000);
    fprintf(stdout, "Azimuth (deg)               [%lf]\n", pf->hdr.azimuth);
    fprintf(stdout, "Zenith Angle (deg)          [%lf]\n", pf->hdr.zenith_ang);
    fprintf(stdout, "Beam FWHM (deg)             [%lf]\n", pf->hdr.beam_FWHM);

    fprintf(stdout, "Length of scan in this file [%lf]\n", pf->hdr.scanlen);
    fprintf(stdout, "Seconds past 00h LST        [%lf]\n", pf->hdr.start_lst);
    fprintf(stdout, "Seconds past 00h UTC        [%lf]\n", pf->hdr.start_sec);

    fprintf(stdout, "Start MJD (whole day)       [%d]\n", pf->hdr.start_day);
    fprintf(stdout, "Scan Number                 [%d]\n", pf->hdr.scan_number);
    fprintf(stdout, "Number of bits per sample   [%d]\n", pf->hdr.nbits);

    fprintf(stdout, "Number of Channels          [%d]\n", pf->hdr.nchan);
    fprintf(stdout, "Number of polarisations     [%d]\n", pf->hdr.npol);
    fprintf(stdout, "Number of spectra per row   [%d]\n", pf->hdr.nsblk);

    fprintf(stdout, "Summed Polarisations? [1/0] [%d]\n", pf->hdr.summed_polns);
    fprintf(stdout, "Receiver Polarisation       [%d]\n", pf->hdr.rcvr_polns);
    fprintf(stdout, "Offset Subint               [%d]\n", pf->hdr.offset_subint);
    fprintf(stdout, "Dwnsmpl fact in time        [%d]\n", pf->hdr.ds_time_fact);
    fprintf(stdout, "Dwnsmpl fact in freq        [%d]\n", pf->hdr.ds_freq_fact);
    fprintf(stdout, "Only Stokes I?              [%d]\n", pf->hdr.onlyI);

    fprintf(stdout, "\nPSRFITS SUBINT:\n");
    fprintf(stdout, "Length of subint (sec)      [%lf]\n", pf->sub.tsubint);
    fprintf(stdout, "Offset (sec)                [%lf]\n", pf->sub.offs);
    fprintf(stdout, "LST (sec)                   [%lf]\n", pf->sub.lst);
    fprintf(stdout, "RA (J2000) (deg)            [%lf]\n", pf->sub.ra);
    fprintf(stdout, "Dec (J2000) (deg)           [%lf]\n", pf->sub.dec);
    fprintf(stdout, "Gal. long. (deg)            [%lf]\n", pf->sub.glon);
    fprintf(stdout, "Gal. lat. (deg)             [%lf]\n", pf->sub.glat);
    fprintf(stdout, "Feed angle (deg)            [%lf]\n", pf->sub.feed_ang);
    fprintf(stdout, "Pos angle of feed (deg)     [%lf]\n", pf->sub.pos_ang);
    fprintf(stdout, "Parallactic angle           [%lf]\n", pf->sub.par_ang);
    fprintf(stdout, "Telescope azimuth           [%lf]\n", pf->sub.tel_az);
    fprintf(stdout, "Telescope zenith angle      [%lf]\n", pf->sub.tel_zen);
    fprintf(stdout, "Bytes per row of raw data   [%d]\n", pf->sub.bytes_per_subint);
    fprintf(stdout, "FITS data typecode          [%d]\n", pf->sub.FITS_typecode);

}

float *create_data_buffer_psrfits( size_t size )
{
    float *ptr = (float *)malloc( size * sizeof(float) );
    return ptr;
}

void init_mpi_psrfits(
        mpi_psrfits *mpf,
        MetafitsMetadata *obs_metadata,
        VoltageMetadata *vcs_metadata,
        int world_size,
        int world_rank,
        int max_sec_per_file,
        int nstokes,
        struct beam_geom *bg,
        char *outfile,
        bool is_writer,
        bool is_coherent )
{
    mpf->ncoarse_chans = world_size;
    mpf->is_writer = is_writer;
    int coarse_chan_idx = vcs_metadata->provided_coarse_chan_indices[0];
    int first_coarse_chan_idx = coarse_chan_idx - world_rank;

    // Populate the PSRFITS header struct for the combined (spliced) output file
    if (is_writer)
        populate_spliced_psrfits_header( &(mpf->spliced_pf), obs_metadata, vcs_metadata,
                first_coarse_chan_idx, mpf->ncoarse_chans, max_sec_per_file, nstokes,
                bg, outfile, is_coherent );

    // Populate the PSRFITS header struct for a single channel
    populate_psrfits_header( &(mpf->coarse_chan_pf), obs_metadata, vcs_metadata, coarse_chan_idx, max_sec_per_file,
            nstokes, bg, outfile, is_coherent );

    // Create MPI vector types designed to splice the coarse channels together
    // correctly during MPI_Gather
    MPI_Type_contiguous( mpf->coarse_chan_pf.hdr.nchan*nstokes, MPI_BYTE, &(mpf->coarse_chan_spectrum) );
    MPI_Type_commit( &(mpf->coarse_chan_spectrum) );

    MPI_Type_vector( mpf->coarse_chan_pf.hdr.nsblk, 1, mpf->ncoarse_chans,
            mpf->coarse_chan_spectrum, &(mpf->total_spectrum_type) );
    MPI_Type_commit( &(mpf->total_spectrum_type) );

    MPI_Type_create_resized( mpf->total_spectrum_type, 0, mpf->coarse_chan_pf.hdr.nchan*nstokes, &(mpf->spliced_type) );
    MPI_Type_commit( &(mpf->spliced_type) );
}

void free_mpi_psrfits( mpi_psrfits *mpf )
{
    MPI_Type_free( &(mpf->spliced_type) );
    MPI_Type_free( &(mpf->total_spectrum_type) );
    MPI_Type_free( &(mpf->coarse_chan_spectrum) );

    free_psrfits( &(mpf->coarse_chan_pf) );

    if (mpf->is_writer)
    {
        free_psrfits( &(mpf->spliced_pf) );
    }
}

void gather_splice_psrfits( mpi_psrfits *mpf, int writer_proc_id )
{
    MPI_Igather(
            mpf->coarse_chan_pf.sub.data, mpf->coarse_chan_pf.hdr.nsblk, mpf->coarse_chan_spectrum,
            mpf->spliced_pf.sub.data, 1, mpf->spliced_type,
            writer_proc_id, MPI_COMM_WORLD, &(mpf->request_data) );

    MPI_Igather(
            mpf->coarse_chan_pf.sub.dat_offsets, mpf->coarse_chan_pf.hdr.nchan, MPI_BYTE,
            mpf->spliced_pf.sub.dat_offsets, mpf->coarse_chan_pf.hdr.nchan, MPI_BYTE,
            writer_proc_id, MPI_COMM_WORLD, &(mpf->request_offsets) );

    MPI_Igather(
            mpf->coarse_chan_pf.sub.dat_scales, mpf->coarse_chan_pf.hdr.nchan, MPI_BYTE,
            mpf->spliced_pf.sub.dat_scales, mpf->coarse_chan_pf.hdr.nchan, MPI_BYTE,
            writer_proc_id, MPI_COMM_WORLD, &(mpf->request_scales) );
}

void wait_splice_psrfits( mpi_psrfits *mpf )
{
    MPI_Wait( &(mpf->request_data),    MPI_STATUS_IGNORE );
    MPI_Wait( &(mpf->request_offsets), MPI_STATUS_IGNORE );
    MPI_Wait( &(mpf->request_scales),  MPI_STATUS_IGNORE );
}
