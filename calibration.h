#ifndef __CALIBRATION_H__
#define __CALIBRATION_H__

#include <inttypes.h>
#include <mwalib.h>
#include <cuComplex.h>

// Calibration solution types
#define NO_CALIBRATION  0
#define RTS             1
#define RTS_BANDPASS    2
#define OFFRINGA        3

#define NDBL_PER_JONES  8

#define CAL_PATHSIZE    4096

struct calibration {
    char  *filename;           // The file that houses the calibration solution
    char  *bandpass_filename;  // The file that houses the RTS bandpass information
    int    nchan;              // The number of channels in the RTS bandpass solutions
    int    cal_type;           // Either RTS or OFFRINGA
    int    offr_chan_num;      // The channel number in the Offringa calibration solution file
    int    ref_ant;            // Reference antenna for calibration phases
    int    cross_terms;        // Include XY and YX of calibration Jones matrices
    double phase_offset;       // Rotate the phase of Y by m*freq + c, where
    double phase_slope;        //   m = phase_slope (rad/Hz)
                               //   c = phase_offset (rad)
};

void get_rts_solution( cuDoubleComplex **D, MetafitsMetadata *cal_metadata,
        const char *caldir, uintptr_t rec_channel );

int read_dijones_file( double **D, double *amp, uintptr_t nant, char *fname );
int read_bandpass_file( cuDoubleComplex ***Jm, cuDoubleComplex ***Jf,
        int chan_width, int nchan, int nant, char *filename );
int read_offringa_gains_file( cuDoubleComplex **antenna_gain, int nant,
        int coarse_chan, char *gains_file, int *order );

uint32_t get_idx_for_vcs_antenna_in_cal( MetafitsMetadata *cal_metadata, MetafitsMetadata *obs_metadata, uint32_t vcs_ant );

#endif
