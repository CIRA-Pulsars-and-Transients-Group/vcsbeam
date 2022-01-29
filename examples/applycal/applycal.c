#include <vcsbeam.h>

int main( int argc, char *argv[] )
{
    // Initialise a VCSBeam context
    vcsbeam_context *vm = vmInit();

    // Print a title message
    vmPrintTitle( vm, "Apply Calibration" );

    // Load an observation's metafits file
    vmLoadObsMetafits( vm, argv[1] );

    // Load the metafits file for the calibration observation
    vmLoadCalMetafits( vm, argv[2] );

    // Prepare the context for reading in a single second and single coarse channel of data
    char *coarse_chan_str        = argv[3];
    char *gps_second_str         = argv[4];
    int   ncoarse_chans          = 1;
    int   nseconds               = 1;
    int   coarse_chan_idx_offset = 0;
    int   gps_second_offset      = 0;
    char *datadir                = argv[5];
    char *calpath                = argv[6];

    // (This also allocates host memory)
    vmBindObsData( vm,
            coarse_chan_str, ncoarse_chans, coarse_chan_idx_offset,
            gps_second_str, nseconds, gps_second_offset,
            datadir );

    // Allocate device memory for voltage data
    vmMallocVDevice( vm );

    // Allocate device and host memory for the calibration solution
    vmMallocDHost( vm );
    vmMallocDDevice( vm );

    // Allocate device and host memory for the antenna mapping index lists
    vmMallocPQIdxsHost( vm );
    vmMallocPQIdxsDevice( vm );

    // Create a lists of rf_input indexes ordered by antenna number (needed for gpu kernels)
    // and upload them to the gpu
    vmSetPolIdxLists( vm );
    vmPushPolIdxLists( vm );

    // Read in the calibration solution (Offringa format)
    bool use_bandpass = false; // (Ignored for Offringa format solutions)
    char *flags_file = NULL; // (Path to file with TileNames to be flagged, if desired)
    vmBindCalibrationData( vm, calpath, CAL_OFFRINGA, use_bandpass, flags_file );
    vmReadCalibration( vm );

    // Read in the voltage data
    vmReadNextSecond( vm );

    // Free device and host memory
    vmFreePQIdxsDevice( vm );
    vmFreePQIdxsHost( vm );
    vmFreeVDevice( vm );
    vmFreeDHost( vm );
    vmFreeDDevice( vm );

    // Finalise
    destroy_vcsbeam_context( vm );

    return 0;
}

