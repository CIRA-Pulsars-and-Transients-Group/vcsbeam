#ifndef __CALIBRATION_H__
#define __CALIBRATION_H__

#include <inttypes.h>
#include <mwalib.h>
#include <cuComplex.h>

// Calibration solution types
#define CAL_NONE      0
#define CAL_RTS       1
#define CAL_OFFRINGA  2

#define NDBL_PER_JONES  8
#define BANDPASS_ROWS_PER_ANT  8

#define CAL_BUFSIZE    4096

struct calibration
{
    char  *metafits;           // Filename of the metafits file
    char  *caldir;             // Location of calibration data
    int    cal_type;           // Either RTS or OFFRINGA
    int    ref_ant;            // Reference antenna for calibration phases
    int    cross_terms;        // Include XY and YX of calibration Jones matrices
    double phase_offset;       // Rotate the phase of Y by m*freq + c, where
    double phase_slope;        //   m = phase_slope (rad/Hz)
                               //   c = phase_offset (rad)
    bool   apply_xy_correction;
};

cuDoubleComplex *get_rts_solution( MetafitsMetadata *cal_metadata,
        MetafitsMetadata *obs_metadata, const char *caldir, uintptr_t coarse_chan_idx );

void read_dijones_file( double **Dd, double *amp, uintptr_t nant, char *fname );
void read_bandpass_file( cuDoubleComplex ***Jm, cuDoubleComplex ***Jf,
        MetafitsMetadata *cal_metadata, char *filename );
int read_offringa_gains_file( cuDoubleComplex **antenna_gain, int nant,
        int coarse_chan, char *gains_file );

void remove_reference_phase( cuDoubleComplex *J, cuDoubleComplex *Jref );
void zero_XY_and_YX( cuDoubleComplex *J );

uint32_t get_idx_for_vcs_antenna_in_cal( MetafitsMetadata *cal_metadata, MetafitsMetadata *obs_metadata, uint32_t vcs_ant );

void pq_phase_correction( uint32_t gpstime, double *phase_slope_rad_per_hz, double *phase_offset_rad );

#endif
