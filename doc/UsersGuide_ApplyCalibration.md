# User's Guide -- Apply Calibration {#usersguideapplycalibration}

[TOC]

## Apply Calibration example progresm

Below is an example program to demonstrate how to apply a calibration solution to a set of voltages using the VCSBeam library.

### Code

**applycal.c:**
```
#include <vcsbeam.h>

void main()
{
    // Initialise a VCSBeam context
    vcsbeam_context *vm = vmInit();

    // Print a Hello World message
    vmPrintTitle( vm, "Apply Calibration" );

    // Finalise
    destroy_vcsbeam_context( vm );
}
```

### Compiling

(There are currently far too many dependencies. This is intended to be cleaned up!!)
```
g++ -o applycal applycal.c \
    -I/opt/cuda/targets/x86_64-linux/include \
    -L/opt/cuda/targets/x86_64-linux/lib \
    -lvcsbeam -lcudart -lcufft \
    -lmwalib \
```

### Running

```
./applycal
```

#### Output (example)

```
------- VCSBeam (v3.3.34_61fcd3e): Apply Calibration -------
```
