# User's Guide -- Hello World {#usersguidehelloworld}

[TOC]

## Hello World program

Below is a basic example program whose purpose is to demonstrate how to compile and run a program using the VCSBeam library.

### Code

```
#include <vcsbeam.h>

void main()
{
    // Initialise a VCSBeam context
    bool use_mpi = false;
    vcsbeam_context *vm = vmInit( use_mpi );

    // Print a Hello World message
    vmPrintTitle( vm, "Hello World!" );

    // Finalise
    destroy_vcsbeam_context( vm );
}
```

### Compiling

(There are currently far too many dependencies. This is intended to be cleaned up!!)
```
g++ -o helloworld helloworld.c \
    -I/opt/cuda/targets/x86_64-linux/include \
    -lvcsbeam -lm -lcudart -lcufft \
    -L/opt/cuda/targets/x86_64-linux/lib \
    -pthread $(mpicc -showme:compile) $(mpicc -showme:link) \
    -lmwalib -lmwa_hyperbeam -lpal \
    -lcudadevrt -lcudart_static -lrt
```

### Running

```
./helloworld
```

#### Output (example)

```
------- VCSBeam (v3.3.34_61fcd3e): Hello World! -------
```
