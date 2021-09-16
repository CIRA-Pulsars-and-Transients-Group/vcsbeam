/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __GEOMETRIC_DELAY_H__
#define __GEOMETRIC_DELAY_H__

#include <cuComplex.h>
#include <mwalib.h>

#include "vcsbeam.h"

/* In order to communicate easily with the GPU, the "phi" array, which
 * contains the complex-valued geometric delay terms, is implemented
 * as a 1D array, which uses the following indexing macro to access
 * the correct term for a given (p)ointing, (a)ntenna, and (c)hannel.
 */

#define PHI_IDX(p,a,c,na,nc)   ((p) * (nc)*(na)  + \
                                (a) * (nc)       + \
                                (c))

typedef struct geometric_delays_t {
    cuDoubleComplex   *phi;
    cuDoubleComplex   *d_phi;
    uintptr_t          npointings;
    uintptr_t          nant;
    uintptr_t          nchan;
    double            *chan_freqs_hz;
    MetafitsMetadata  *obs_metadata;
} geometric_delays;

/* Calculate the geometric delay (in radians) for the given pointings
 */
void calc_all_geometric_delays(
        geometric_delays  *gdelays,
        struct beam_geom  *beam_geom_vals );

void calc_geometric_delays(
        struct beam_geom  *beam_geom_vals,
        uint32_t           freq_hz,
        MetafitsMetadata  *obs_metadata,
        cuDoubleComplex   *phi );

/* Allocates memory for the geometric delay arrays ("phi") on both host and device.
 * Free with free_geometric_delays()
 */
void create_geometric_delays(
        geometric_delays  *gdelays,
        MetafitsMetadata  *obs_metadata,
        VoltageMetadata   *vcs_metadata,
        uintptr_t          coarse_chan,
        uintptr_t          npointings );


/* Free memory allocated with create_geometric_delays()
 */
void free_geometric_delays( geometric_delays *gdelays );

/* Copy host memory block to device
 */
void push_geometric_delays_to_device( geometric_delays *gdelays );

double calc_array_factor(
        MetafitsMetadata *obs_metadata,
        uint32_t          freq_hz,
        struct beam_geom *bg1,
        struct beam_geom *bg2 );

void calc_beam_geom(
        double            ras_hours,
        double            decs_degs,
        double            mjd,
        struct beam_geom  *bg );

void dec2hms( char *out, double in, int sflag );
void utc2mjd( char *, double *, double * ); // "2000-01-01T00:00:00" --> MJD_int + MJD_fraction
void mjd2lst( double, double * );

double parse_dec( char* ); // "01:23:45.67" --> Dec in degrees
double parse_ra( char* );  // "01:23:45.67" --> RA  in degrees


#endif
