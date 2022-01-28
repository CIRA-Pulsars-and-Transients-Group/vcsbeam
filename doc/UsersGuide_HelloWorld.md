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

```
gcc -o helloworld helloworld.c -lvcsbeam
```

### Running

```
./helloworld
```

#### Output

```
...
```
