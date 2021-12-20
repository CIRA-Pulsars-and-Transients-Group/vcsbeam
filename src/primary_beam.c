/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cuComplex.h>

#include <mwalib.h>
#include <mwa_hyperbeam.h>

#include "vcsbeam.h"

#define NCOMPLEXELEMENTS 4

const double Isky[] = { 0.5, 0. , 0. ,  0. , 0. , 0. ,  0.5, 0. };
const double Qsky[] = { 0.5, 0. , 0. ,  0. , 0. , 0. , -0.5, 0. };
const double Usky[] = { 0. , 0. , 0.5,  0. , 0.5, 0. ,  0. , 0. };
const double Vsky[] = { 0. , 0. , 0. , -0.5, 0  , 0.5,  0. , 0. };
const double *sky[] = { Isky, Qsky, Usky, Vsky };

void vmCalcB(
        vcsbeam_context   *vm,
        beam_geom         *beam_geom_vals )
/* Calculate the required beam matrices for the given pointings. The
 * calculated Jones matrices are stored in pb->B
 *
 * Only those configurations of live (i.e. non-dead) dipoles which exist in
 * the array are calculated, in order to save calculation time.
 */
{
    // Shorthand variable for where to put the answer
    primary_beam *pb = &vm->pb;

    // Calculate some array sizes
    uintptr_t nant      = pb->nant;
    uintptr_t npol      = pb->npol; // = 4 (XX, XY, YX, YY)
    uintptr_t nrf_input = pb->obs_metadata->num_rf_inputs;

    // Set up the output array of beam matrices ("B"; see Eq. (30) in Ord et al. (2019))
    // Layout is B[pointing][antenna][pol] where 0 <= pol < npol=4
    uintptr_t p, ant;

    // Make temporary array that will hold jones matrices for specific configurations
    cuDoubleComplex *configs[NCONFIGS];
    int config_idx;
    for (config_idx = 0; config_idx < NCONFIGS; config_idx++)
        configs[config_idx] = NULL;

    // Normalise to zenith
    int zenith_norm = 1;

    // Prepare a parallactic angle correction rotation matrix
    double P[npol];

    // We will be looping over the antennas/pols and checking each dipole
    // configuration to see if it's been already calculated. However, we
    // still need to know how many antennas/pols (= RF inputs) there are
    uintptr_t rf_input;

    double az, za; // (az)imuth and (z)enith (a)ngle in radians

    // Loop through the pointings and calculate the primary beam matrices
    for (p = 0; p < pb->npointings; p++)
    {
        az = beam_geom_vals[p].az;
        za = PIBY2 - beam_geom_vals[p].el;

        // Calculate the parallactic angle correction for this pointing
        parallactic_angle_correction( P, MWA_LATITUDE_RADIANS, az, za );

        for (rf_input = 0; rf_input < nrf_input; rf_input++)
        {
            // Get the antenna from the rf_input
            ant = pb->obs_metadata->rf_inputs[rf_input].ant;

            // Assume that the delay config for the 'Y' pol matches that of the corresponding 'X'
            if (*(pb->obs_metadata->rf_inputs[rf_input].pol) == 'Y')
                continue;

            // Check this RF input's dipole configuration
            config_idx = hash_dipole_config( pb->amps[rf_input] );
            if (config_idx == NCONFIGS - 1)
                continue; // All dipoles are dead -- skip this tile

            if (config_idx == MANY_DEAD_DIPOLES) // If lots of dipoles are dead, call Hyperbeam
            {
                // Use the 'dead' configuration temporarily
                config_idx = DEAD_CONFIG;
                configs[config_idx] = (cuDoubleComplex *)calc_jones(
                        pb->beam, az, za, pb->freq_hz, pb->delays[rf_input], pb->amps[rf_input], zenith_norm );
            }
            else if (configs[config_idx] == NULL) // Call Hyperbeam if this config hasn't been done yet
            {
                // Get the calculated FEE Beam (using Hyperbeam)
                configs[config_idx] = (cuDoubleComplex *)calc_jones(
                        pb->beam, az, za, pb->freq_hz, pb->delays[rf_input], pb->amps[rf_input], zenith_norm );

                // Apply the parallactic angle correction
#ifdef DEBUG
if (config_idx == 0)
{
    fprintf( stderr, "B       = " ); fprintf_complex_matrix( stderr, configs[config_idx] );
    fprintf( stderr, "P       = " ); fprintf( stderr, "[ %lf, %lf; %lf, %lf ]\n", P[0], P[1], P[2], P[3]  );
}
#endif
                mult2x2d_CxR( configs[config_idx], P, configs[config_idx] );
//fprintf( stderr, "after  pa correction: B = " );
//fprintf_complex_matrix( stderr, configs[config_idx] );
            }

            // Copy the answer into the B matrix (for this antenna)
            cp2x2( configs[config_idx], &(pb->B[PB_IDX(p, ant, 0, nant, npol)]) );
        }
    }
}

void vmCreatePrimaryBeam( vcsbeam_context *vm )
/* Allocates memory for the primary beam matrices ("B")
 * (see Eq. (30) in Ord et al. (2019))
 */
{
    // Calculate some array sizes
    vm->pb.npointings = vm->npointing;
    vm->pb.nant = vm->obs_metadata->num_ants;
    vm->pb.npol = vm->obs_metadata->num_visibility_pols; // = 4 (XX, XY, YX, YY)

    size_t size = vm->pb.npointings * vm->pb.nant * vm->pb.npol * sizeof(cuDoubleComplex);

    // Allocate memory
    vm->pb.B = (cuDoubleComplex *)malloc( size );

    vm->pb.beam = NULL;
    vm->pb.beam = new_fee_beam( HYPERBEAM_HDF5 );

    create_delays_amps_from_metafits( vm->obs_metadata, &(vm->pb.delays), &(vm->pb.amps) );

    vm->pb.freq_hz = vm->obs_metadata->metafits_coarse_chans[vm->coarse_chan_idx].chan_centre_hz;

    vm->pb.obs_metadata = vm->obs_metadata;
}

void free_primary_beam( primary_beam *pb )
/* Frees the memory allocated in malloc_primary_beam()
 */
{
    free( pb->B );
    free_delays_amps( pb->obs_metadata, pb->delays, pb->amps );
    free_fee_beam( pb->beam );
}


void create_delays_amps_from_metafits( MetafitsMetadata *obs_metadata, uint32_t ***delays, double ***amps )
/* Costruct "delay" and "amp" arrays required by Hyperbeam.
 * Delays can be brought directly across from the MetafitsMetadata struct.
 * The "conversion" from delays to amps is straightforward: All amps are "1.0"
 * unless the delay is 32, in which case the amp should be "0.0".
 * (The value 32 is the number assigned to a faulty or broken dipole.)
 *
 * Both the delays array (input) and the amps array (output) should
 * have dimensions [ninputs][ndipoles]
 *
 * This function allocates memory for both delays and amps. This can be freed
 * by a call to free_delays_amps().
 */
{
    int ndipoles = obs_metadata->num_delays;
    int ninputs  = obs_metadata->num_rf_inputs;

    *delays = (uint32_t **)malloc( ninputs * sizeof(uint32_t *) );
    *amps = (double **)malloc( ninputs * sizeof(double *) );

    int input, dipole;

    // For the delays, we will just point each delay array to the
    // corresponding array in the metafits, but for amps, we need
    // to allocate new memory for each set of ndipoles (=16) numbers
    // and populate them
    for (input = 0; input < ninputs; input++)
    {
        (*delays)[input] = obs_metadata->rf_inputs[input].dipole_delays;

        (*amps)[input] = (double *)malloc( ndipoles * sizeof(double) );
        for (dipole = 0; dipole < ndipoles; dipole++)
        {
            (*amps)[input][dipole] = ((*delays)[input][dipole] == 32 ? 0.0 : 1.0);
        }
    }
}


void free_delays_amps( MetafitsMetadata *obs_metadata, uint32_t **delays, double **amps )
/* Frees the memory allocated in create_delays_amps_from_metafits().
 */
{
    int ninputs  = obs_metadata->num_rf_inputs;
    int input;
    for (input = 0; input < ninputs; input++)
        free( amps[input] );

    free( amps );
    free( delays );
}


int hash_dipole_config( double *amps )
/* In order to avoid recalculating the FEE beam for repeated dipole
 * configurations, we have to keep track of which configurations have already
 * been calculated. We do this through a boolean array, and this function
 * converts dipole configurations into indices of this array. In other words,
 * this function _assigns meaning_ to the array.
 *
 * Since dead dipoles are relatively rare, we only consider configurations
 * in which up to two dipoles are dead. Any more than that and the we can
 * recalculate the Jones matrix with minimal entropy. In this case, this
 * function returns MANY_DEAD_DIPOLES. The other remaining cases are:
 *
 *   0  dead dipoles = 1   configuration
 *   1  dead dipole  = 16  configurations
 *   2  dead dipoles = 120 configurations
 *   16 dead dipoles = 1   configuration
 *
 * for a grand total of 138 indices. They are ordered as follows:
 *
 *  idx  configuration (0=dead, 1=live)
 *   0   [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   1   [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   2   [1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   3   [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   ...
 *   16  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
 *   17  [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   18  [0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   19  [0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   ...
 *   31  [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
 *   32  [1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   32  [1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   ...
 *   136 [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0]
 *   136 [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
 *   137 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 *
 *   This function defines "dead" as amps=0.0, "live" otherwise.
 */
{
    // The return value = the "hashed" array index
    int idx;
    int d1 = 0, d2 = 0; // The "locations" of the (up to) two dead dipoles

    // Count the number of dead dipoles
    int ndead = 0;
    int ndipoles = 16;
    int i;
    for (i = 0; i < ndipoles; i++)
        ndead += (amps[i] == 0.0);

    // Convert configuration into an array index (the "hash" function)
    switch (ndead)
    {
        case 0:
            idx = 0;
            break;
        case 1:
            for (i = 0; i < ndipoles; i++)
            {
                if (amps[i] == 0.0)
                {
                    d1 = i;
                    break;
                }
            }
            idx = d1 + 1; // The "1-dead-dipole" configs start at idx 1
            break;
        case 2:
            // Find the locations of the two dead dipoles
            d1 = -1;
            for (i = 0; i < ndipoles; i++)
            {
                if (amps[i] == 0.0)
                {
                    if (d1 == -1)
                        d1 = i;
                    else
                    {
                        d2 = i;
                        break;
                    }
                }
            }
            // The hash function for two dead dipoles
            // (The "2-dead-dipole" configs start at idx 17
            idx = 16*d1 - ((d1 + 2)*(d1 + 1))/2 + d2 + 17;
            break;
        case 16:
            idx = 137;
            break;
        default: // any configuration with >2 dead dipoles
            idx = MANY_DEAD_DIPOLES;
            break;
    }

    return idx;
}


void parallactic_angle_correction(
    double *Ppa,  // output rotation matrix
    double lat,   // observing latitude (radians)
    double az,    // azimuth angle (radians)
    double za)    // zenith angle (radians)
{
    /* This parallactic angle correction is intended to be applied to the
     * primary beam matrices delivered by the FEE beam code
     * The FEE beam Jones matrix (according to Sokolowski et al. 2017) is
     *   [q] = [ qθ  qφ ] [θ]
     *   [p]   [ pθ  pφ ] [φ]
     * However, we would like a matrix that goes from (y,x) to (q,p):
     *   [q] = [ qy  qx ] [y]
     *   [p]   [ py  px ] [x]
     * This means we will have to multiply the FEE matrix on the right by
     *   [q] = [ qθ  qφ ] [ θy θx ] [y]
     *   [p]   [ pθ  pφ ] [ φy φx ] [x]
     * The second matrix on the RHS above is the parallactic angle correction.
     * For these specific coordinate systems,
     *   [ θy θx ] = [ -sin(χ)  -cos(χ) ]  (see docs for sign conventions used here)
     *   [ φy φx ]   [ -cos(χ)   sin(χ) ]
     * Thus, the parallactic angle correction computed here is the transformation
     * from (y,x) to (θ,φ)
     */

    double el = PIBY2 - za;

    double ha, dec;
    palDh2e(az, el, lat, &ha, &dec);
    double chi = palPa( ha, dec, lat );

    double schi = sin(chi);
    double cchi = cos(chi);

    Ppa[0] = -cchi;
    Ppa[1] = -schi;
    Ppa[2] =  schi;
    Ppa[3] = -cchi;
}

void calc_normalised_beam_response( FEEBeam *beam, double az, double za, double freq_hz, uint32_t *delays, double *amps, double *IQUV, cuDoubleComplex **J, bool apply_pa_correction )
{
    cuDoubleComplex JH[NCOMPLEXELEMENTS];
    cuDoubleComplex sky_x_JH[NCOMPLEXELEMENTS];
    cuDoubleComplex coherency[NCOMPLEXELEMENTS];

    // Calculate the primary beam for this channel, in this direction
    int zenith_norm = 1;
    *J = (cuDoubleComplex *)calc_jones( beam, az, za, freq_hz, delays, amps, zenith_norm );

    // Optionally apply the parallactic angle correction
    double P[NCOMPLEXELEMENTS]; // (Real-valued) parallactic angle correction matrix
    if (apply_pa_correction)
    {
        // Calculate the parallactic angle correction
        parallactic_angle_correction( P, MWA_LATITUDE_RADIANS, az, za );

        mult2x2d_RxC( P, *J, *J );
    }

    // Convert this jones matrix to Stokes parameters for I, Q, U, V skies
    int stokes;
    for (stokes = 0; stokes < 4; stokes++)
    {
        calc_hermitian( *J, JH );
        mult2x2d( (cuDoubleComplex *)sky[stokes], JH, sky_x_JH );
        mult2x2d( *J, sky_x_JH, coherency );

        if (stokes == 0) // Stokes I
            IQUV[0] = 0.5*cuCreal( cuCadd( coherency[0], coherency[3] ) );
        else if (stokes == 1) // Stokes Q
            IQUV[1] = 0.5*cuCreal( cuCsub( coherency[0], coherency[3] ) );
        else if (stokes == 2) // Stokes U
            IQUV[2] =  cuCreal( coherency[1] );
        else  // if (stokes == 3) // Stokes V
            IQUV[3] = -cuCimag( coherency[1] );
    }
}

