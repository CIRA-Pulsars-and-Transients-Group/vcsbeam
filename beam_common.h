/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef BEAM_COMMON_H
#define BEAM_COMMON_H

#include <inttypes.h>
#include "mwa_hyperbeam.h"
#include "calibration.h"
#include <cuComplex.h>
#include <mwalib.h>

#define MAX_POLS   4

#define N_COPOL    2
#define R2C_SIGN   -1.0
#define NDELAYS    16

#define MWA_LAT -26.703319        /* Array latitude. degrees North */
#define MWA_LON 116.67081         /* Array longitude. degrees East */
#define MWA_HGT 377               /* Array altitude. meters above sea level */

// A structure to read in all the relevant info from the observation metafits
// file.
struct metafits_info {
    double      tile_pointing_ra;
    double      tile_pointing_dec;
    double      tile_pointing_az;
    double      tile_pointing_el;
    float      *N_array;
    float      *E_array;
    float      *H_array;
    float      *cable_array;
    int        *flag_array;
    double     *weights_array;
    short int  *antenna_num;
    char      **tilenames;
    int         ninput;
    int         chan_width;
    int       **delays;
    double    **amps;
    char       *date_obs;
    int         exposure;
};

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

struct make_beam_opts {
    // Variables for required options
    unsigned long int  begin;         // GPS time -- when to start beamforming
    unsigned long int  end;           // GPS time -- when to stop beamforming
    char              *pointings;     // pointing list"dd:mm:ss_hh:mm:ss, dd:mm:ss_hh:mm:ss"
    char              *datadir;       // The path to where the recombined data live
    char              *metafits;      // filename of the metafits file
    char              *cal_metafits;  // filename of the metafits file
    int                rec_channel;   // 0 - 255 receiver 1.28MHz channel
    long int           frequency;     // = rec_channel expressed in Hz

    // Variables for MWA/VCS configuration
    char              *custom_flags;  // Use custom list for flagging antennas

    // Output options
    int                out_incoh;     // Default = PSRFITS (incoherent) output turned OFF
    int                out_coh;       // Default = PSRFITS (coherent)   output turned OFF
    int                out_vdif;      // Default = VDIF                 output turned OFF
    int                out_uvdif;     // Default = upsampled VDIF       output turned OFF
    int                out_bf;        // Default = beamform over all (non-flagged) antennas
    int                out_ant;       // The antenna number (0-127) to write out if out_bf = 0

    // Other options
    char              *synth_filter;  // Which synthesis filter to use
    int                out_summed;    // Default = output only Stokes I output turned OFF
    int                max_sec_per_file;    // Number of seconds per fits files
    float              gpu_mem  ;     // Default = -1.0. If -1.0 use all GPU mem
};

/* Running get_jones from within make_beam */
void get_jones(
        // an array of pointings [pointing][ra/dec][characters]
        int                    npointing, // number of pointings
        VoltageMetadata*       volt_metadata,
        MetafitsMetadata      *metafits_metadata,
        MetafitsMetadata*      cal_metadata,
        int                    coarse_chan_idx,
        struct                 calibration *cal,
        cuDoubleComplex      **M,
        cuDoubleComplex     ***Jf,
        float                  samples_per_sec,
        FEEBeam               *beam,
        uint32_t             **delays,
        double               **amps,
        struct beam_geom       beam_geom_vals[],
        cuDoubleComplex      ****complex_weights_array,  // output
        cuDoubleComplex      ****invJi                   // output
);

void calc_beam_geom(
        char              pointing_array[][2][64],
        int               npointing,
        double            mjd,
        struct beam_geom  bg[] );

void create_delays_amps_from_metafits( MetafitsMetadata *metafits_metadata, uint32_t ***delays, double ***amps );
void free_delays_amps( MetafitsMetadata *metafits_metadata, uint32_t **delays, double **amps );

void create_antenna_lists( MetafitsMetadata *metafits_metadata, uint32_t *polX_idxs, uint32_t *polY_idxs );
void remove_reference_phase( cuDoubleComplex **M, int ref_ant, int nant );
void zero_XY_and_YX( cuDoubleComplex **M, int nant );

int calcEjones_analytic(cuDoubleComplex response[MAX_POLS], // pointer to 4-element (2x2) voltage gain Jones matrix
               const long freq, // observing freq (Hz)
               const float lat, // observing latitude (radians)
               const float az0, // azimuth & zenith angle of tile pointing
               const float za0,
               const float az, // azimuth & zenith angle to sample
               const float za);


void parallactic_angle_correction_analytic(
    double *P,    // output rotation matrix
    double lat,   // observing latitude (radians)
    double az,    // azimuth angle (radians)
    double za);   // zenith angle (radians)

void parallactic_angle_correction_fee2016(
    double *P,    // output rotation matrix
    double lat,   // observing latitude (radians)
    double az,    // azimuth angle (radians)
    double za);   // zenith angle (radians)

int hash_dipole_config( double * );

void get_metafits_info( char *metafits, struct metafits_info *mi, unsigned int chan_width );
void destroy_metafits_info( struct metafits_info *mi );

void int8_to_uint8(int n, int shift, char * to_convert);
void float2int8_trunc(float *f, int n, float min, float max, int8_t *i);
void float_to_unit8(float * in, int n, int8_t *out);

void flatten_bandpass(
        int nstep,
        int nchan,
        int npol,
        void *data);

void read_data( char *filename, uint8_t *data, int nbytes );

void dec2hms( char *out, double in, int sflag );
void utc2mjd( char *, double *, double * ); // "2000-01-01T00:00:00" --> MJD_int + MJD_fraction
void mjd2lst( double, double * );

double parse_dec( char* ); // "01:23:45.67" --> Dec in degrees
double parse_ra( char* );  // "01:23:45.67" --> RA  in degrees

/**** MATRIX OPERATIONS ****/

void cp2x2(cuDoubleComplex *Min, cuDoubleComplex *Mout);
void inv2x2(cuDoubleComplex *Min, cuDoubleComplex *Mout);
void inv2x2d(double *Min, double *Mout);
void inv2x2S(cuDoubleComplex *Min, cuDoubleComplex **Mout);
void mult2x2d(cuDoubleComplex *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout);
void mult2x2d_RxC(double *M1, cuDoubleComplex *M2, cuDoubleComplex *Mout);
void conj2x2(cuDoubleComplex *M, cuDoubleComplex *Mout);
double norm2x2(cuDoubleComplex *M, cuDoubleComplex *Mout);

#endif
