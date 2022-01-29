#include <vcsbeam.h>

int main()
{
    // Initialise a VCSBeam context
    vcsbeam_context *vm = vmInit();

    // Print a Hello World message
    vmPrintTitle( vm, "Apply Calibration" );

    // Finalise
    destroy_vcsbeam_context( vm );

    return 0;
}

