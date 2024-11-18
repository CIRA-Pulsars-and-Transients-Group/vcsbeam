/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <mwalib.h>
#include <mwa_hyperbeam.h>

#include "vcsbeam.h"
#include "gpu_macros.h"

#define NCOMPLEXELEMENTS 4

/**
 * The coherency matrix for a pure Stokes I sky.
 */
const double Isky[] = { 0.5, 0. , 0. ,  0. , 0. , 0. ,  0.5, 0. };
/**
 * The coherency matrix for a pure Stokes Q sky.
 */
const double Qsky[] = { 0.5, 0. , 0. ,  0. , 0. , 0. , -0.5, 0. };
/**
 * The coherency matrix for a pure Stokes U sky.
 */
const double Usky[] = { 0. , 0. , 0.5,  0. , 0.5, 0. ,  0. , 0. };
/**
 * The coherency matrix for a pure Stokes V sky.
 */
const double Vsky[] = { 0. , 0. , 0. , -0.5, 0  , 0.5,  0. , 0. };
/**
 * The set of coherency matrices for pure Stokes I, Q, U, and V skies.
 */
const double *sky[] = { Isky, Qsky, Usky, Vsky };

/**
 * Calculates the beam model Jones matrices for the given pointings.
 *
 * @param vm The VCSBeam context struct
 * @param beam_geom_vals An array of beam_geom objects containing the pointing
 *                       information
 *
 * This function uses Hyperbeam to calculate the beam model Jones matrices
 * for each pointing in `beam_geom_vals`.
 *
 * Only those configurations of live (i.e. non-dead) dipoles which exist in
 * the array are calculated, in order to save calculation time.
 *
 * The parallactic angle correction is applied, so that the final matrices,
 * which are stored in `pb&rarr;B`, are in the basis
 * \f[
 *     {\bf B} =
 *         \begin{bmatrix} B_{qx} & B_{qy} \\ B_{px} & B_{py} \end{bmatrix}.
 * \f]
 */




void handle_hyperbeam_error(char file[], int line_num, const char function_name[]) {
    int err_length = hb_last_error_length(); 
    char *err = malloc(err_length * sizeof(char));
    int err_status = hb_last_error_message(err, err_length);
    if (err_status == -1) {
        printf("Something really bad happened!\n");
        exit(EXIT_FAILURE);
    }
    printf("File %s:%d: hyperbeam error in %s: %s\n", file, line_num, function_name, err);

    exit(EXIT_FAILURE);
}



void vmCalcB(
        vcsbeam_context   *vm,
        beam_geom         *beam_geom_vals )
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
    gpuDoubleComplex *configs[NCONFIGS];
    int config_idx;
    for (config_idx = 0; config_idx < NCONFIGS; config_idx++)
        configs[config_idx] = NULL;

    double * tempJones;
    tempJones = malloc(8*sizeof(double));

    // Normalise to zenith
    int zenith_norm = 1;

    // Prepare a parallactic angle correction rotation matrix
    double P[npol];

    // We will be looping over the antennas/pols and checking each dipole
    // configuration to see if it's been already calculated. However, we
    // still need to know how many antennas/pols (= RF inputs) there are
    uintptr_t rf_input;

    double az, za; // (az)imuth and (z)enith (a)ngle in radians

    uint8_t iauOrder=0; //Boolean 0= don't set Jones matrix to be iau order but instead mwa order in CalcJones

    uint32_t numAmps=16; //number of dipole gains used (16 or 32)

    int32_t errInt; //exit code integer for calcJones

    double * arrayLatitudeRad=NULL;

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
                errInt = fee_calc_jones(pb->beam, az, za, pb->freq_hz, pb->delays[rf_input], pb->amps[rf_input],numAmps, zenith_norm, arrayLatitudeRad ,iauOrder , tempJones);
                configs[config_idx]=(gpuDoubleComplex *)(tempJones);
            }
            else if (configs[config_idx] == NULL) // Call Hyperbeam if this config hasn't been done yet
            {
                // Get the calculated FEE Beam (using Hyperbeam)
                errInt = fee_calc_jones(pb->beam, az, za, pb->freq_hz, pb->delays[rf_input], pb->amps[rf_input],numAmps, zenith_norm, arrayLatitudeRad ,iauOrder , tempJones);
                configs[config_idx]=(gpuDoubleComplex *)(tempJones);

                // Apply the parallactic angle correction
#ifdef DEBUG
                if (config_idx == 0)
                {
                    fprintf( stderr, "Bhb     = " ); fprintf_complex_matrix( stderr, configs[config_idx] );
                    fprintf( stderr, "Ppa     = " ); fprintf( stderr, "[ %lf, %lf; %lf, %lf ]\n", P[0], P[1], P[2], P[3]  );
                }
#endif
                mult2x2d_CxR( configs[config_idx], P, configs[config_idx] );
#ifdef DEBUG
                if (config_idx == 0)
                {
                    fprintf( stderr, "B       = " ); fprintf_complex_matrix( stderr, configs[config_idx] );
                }
#endif
            }

            if (errInt !=0)
            {             
                handle_hyperbeam_error("Primary Beam",183, "fee_calc_jones");   
                exit(EXIT_FAILURE);
            }

            // Copy the answer into the B matrix (for this antenna)
            cp2x2( configs[config_idx], &(pb->B[PB_IDX(p, ant, 0, nant, npol)]) );
        }
    }
    free(tempJones);
}

/**
 * Allocates memory for the beam model Jones matrices.
 *
 * The resulting array (`vm&rarr;B`) has dimensions
 * \f$N_b \times N_a \times N_p \times N_p\f$.
 */
void vmCreatePrimaryBeam( vcsbeam_context *vm )
{
    int32_t  errInt; //new_fee_beam Error integer
    
    // Calculate some array sizes
    vm->pb.npointings = vm->npointing;
    vm->pb.nant = vm->obs_metadata->num_ants;
    vm->pb.npol = vm->obs_metadata->num_visibility_pols; // = 4 (XX, XY, YX, YY)

    size_t size = vm->pb.npointings * vm->pb.nant * vm->pb.npol * sizeof(gpuDoubleComplex);

    // Allocate memory
    vm->pb.B = (gpuDoubleComplex *)malloc( size );

    vm->pb.beam = NULL;
    errInt= new_fee_beam( HYPERBEAM_HDF5, &vm->pb.beam );

    if (errInt != 0)
    {
        handle_hyperbeam_error("Primary Beam", 213, "new_fee_beam");
        exit(EXIT_FAILURE);
    }

    create_delays_amps_from_metafits( vm->obs_metadata, &(vm->pb.delays), &(vm->pb.amps) );

    vm->pb.freq_hz = vm->obs_metadata->metafits_coarse_chans[vm->coarse_chan_idx].chan_centre_hz;

    vm->pb.obs_metadata = vm->obs_metadata;
}

/**
 * Frees the memory allocated in vmCreatePrimaryBeam()
 *
 * @todo Change name to vm...
 */
void free_primary_beam( primary_beam *pb )
{
    free( pb->B );
    free_delays_amps( pb->obs_metadata, pb->delays, pb->amps );
    free_fee_beam( pb->beam );
}


/**
 * Constructs "delay" and "amp" arrays required by Hyperbeam.
 *
 * @param[in]  obs_metadata A MetafitsMetadata object containing the delays
 *                          and amps.
 * @param[out] delays       A pointer to an array of delay values, with layout
 *                          \f$N_i \times N_d\f$
 * @param[out] amps         A pointer to an array of amp values, with layout
 *                          \f$N_i \times N_d\f$
 *
 * Delays can be brought directly across from the MetafitsMetadata struct.
 * The "conversion" from delays to amps is straightforward: All amps are "1.0"
 * unless the delay is 32, in which case the amp should be "0.0".
 * (The value 32 is the number assigned to a faulty or broken dipole.)
 *
 * This function allocates memory for both delays and amps. This can be freed
 * by a call to free_delays_amps().
 */
void create_delays_amps_from_metafits( MetafitsMetadata *obs_metadata, uint32_t ***delays, double ***amps )
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


/**
 * Frees the memory allocated in create_delays_amps_from_metafits().
 *
 * @param obs_metadata A MetafitsMetadata object
 * @param delays       The array of delay values to be freed
 * @param amps         The array of amp values to be freed
 *
 * Here, the MetafitsMetadata is used to find \f$N_i\f$ (one of the dimensions
 * of `delays` and `amps`).
 */
void free_delays_amps( MetafitsMetadata *obs_metadata, uint32_t **delays, double **amps )
{
    int ninputs  = obs_metadata->num_rf_inputs;
    int input;
    for (input = 0; input < ninputs; input++)
        free( amps[input] );

    free( amps );
    free( delays );
}


/**
 * Implements a hash lookup for dead/live dipole configurations.
 *
 * @param amps The array of amps whose configuration is to be hashed
 * @return The hashed index for this configuration
 *
 * In order to avoid recalculating the FEE beam for repeated dipole
 * configurations, we have to keep track of which configurations have already
 * been calculated. We do this through a boolean array, and this function
 * converts dipole configurations into indices of this array. In other words,
 * this function _assigns meaning_ to the array.
 *
 * Since dead dipoles are relatively rare, we only consider configurations
 * in which up to two dipoles are dead. Any more than that and the we can
 * recalculate the Jones matrix with minimal entropy. In this case, this
 * function returns `MANY_DEAD_DIPOLES`. The other remaining cases are:
 *
 *  - 0  dead dipoles = 1   configuration
 *  - 1  dead dipole  = 16  configurations
 *  - 2  dead dipoles = 120 configurations
 *  - 16 dead dipoles = 1   configuration
 *
 * for a grand total of 138 indices. They are ordered as follows:
 * ```
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
 *   33  [1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 *   ...
 *   135 [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0]
 *   136 [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
 *   137 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 * ```
 *
 *  This function defines "dead" as amps=0.0, "live" otherwise.
 */
int hash_dipole_config( double *amps )
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


/**
 * Calculates the parallactic angle correction matrix,
 * \f${\bf P}_\text{pa}\f$.
 *
 * @param[out] Ppa The output rotation matrix
 * @param[in]  lat The observing latitude (radians)
 * @param[in]  az  The azimuth angle (radians)
 * @param[in]  za  The zenith angle (radians)
 *
 * This function computes the parallactic angle correction matrix in the
 * basis (see [Parallactic angle correction](@ref parallacticangle)):
 * \f[
 *     {\bf P}_\text{pa} =
 *     \begin{bmatrix}
 *         P_{\theta x} & P_{\theta y} \\
 P_{\phi x}   & P_{\phi y}
 \end{bmatrix}.
 * \f]
 */
void parallactic_angle_correction(
        double *Ppa,  // output rotation matrix
        double lat,   // observing latitude (radians)
        double az,    // azimuth angle (radians)
        double za)    // zenith angle (radians)
{

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

/**
 * Calculates the normalised primary beam response to a Stokes I, Q, U, or
 * V sky.
 *
 * @param      beam A Hyperbeam object
 * @param[in]  az The azimuth (radians)
 * @param[in]  za The zenith angle (radians)
 * @param[in]  freq_hz The frequency (Hz)
 * @param[in]  delays The dipole delays that define a primary beam pointing
 * @param[in]  amps The dipole amplitudes that define the live/dead
 *             configuration
 * @param[out] IQUV The beam response
 * @param[out] J A pointer to the beam model Jones matrix
 * @param      apply_pa_correction Whether to apply the parallactic angle
 *             correction to `J`
 *
 * `IQUV` must point to already-allocated memory with size \f$N_s\f$.
 *
 * This function will allocate memory for `J`, which must be freed by the
 * caller.
 */
void calc_normalised_beam_response( FEEBeam *beam, double az, double za, double freq_hz, uint32_t *delays, double *amps, double *IQUV, gpuDoubleComplex *J, bool apply_pa_correction )
{

    uint8_t iauOrder=0; //Boolean 0= don't set Jones matrix to be iau order but instead mwa order in CalcJones

    uint32_t numAmps=16; //number of dipole gains used (16 or 32)

    int32_t errInt=0; //exit code integer for calcJones

    double * arrayLatitudeRad;             
    arrayLatitudeRad=NULL;

    gpuDoubleComplex JH[NCOMPLEXELEMENTS];
    gpuDoubleComplex sky_x_JH[NCOMPLEXELEMENTS];
    gpuDoubleComplex coherency[NCOMPLEXELEMENTS];

    // Calculate the primary beam for this channel, in this direction
    int zenith_norm = 1;
    errInt = fee_calc_jones( beam, az, za, freq_hz, delays, amps,numAmps, zenith_norm, arrayLatitudeRad ,iauOrder , (double *)J );

    if (errInt !=0)
    {
        handle_hyperbeam_error("primary beam", __LINE__ ,"fee_calc_jones");
        exit(EXIT_FAILURE);
    }

    // Optionally apply the parallactic angle correction
    double P[NCOMPLEXELEMENTS]; // (Real-valued) parallactic angle correction matrix
    if (apply_pa_correction)
    {
        // Calculate the parallactic angle correction
        parallactic_angle_correction( P, MWA_LATITUDE_RADIANS, az, za );

        mult2x2d_RxC( P, J, J );
    }

    // Convert this jones matrix to Stokes parameters for I, Q, U, V skies
    int stokes;
    for (stokes = 0; stokes < 4; stokes++)
    {
        calc_hermitian( J, JH );
        mult2x2d( (gpuDoubleComplex *)sky[stokes], JH, sky_x_JH );
        mult2x2d( J, sky_x_JH, coherency );

        if (stokes == 0) // Stokes I
            IQUV[0] = 0.5*gpuCreal( gpuCadd( coherency[0], coherency[3] ) );
        else if (stokes == 1) // Stokes Q
            IQUV[1] = 0.5*gpuCreal( gpuCsub( coherency[0], coherency[3] ) );
        else if (stokes == 2) // Stokes U
            IQUV[2] =  gpuCreal( coherency[1] );
        else  // if (stokes == 3) // Stokes V
            IQUV[3] = -gpuCimag( coherency[1] );
    }
}

