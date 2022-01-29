#include <vcsbeam.h>

int main()
{
    // Initialise a VCSBeam context
    vcsbeam_context *vm = vmInit();

    // Print a Hello World message
    vmPrintTitle( vm, "Hello World!" );

    // Finalise
    destroy_vcsbeam_context( vm );

    return 0;
}

