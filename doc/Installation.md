# Installation {#installationguide}

VCSBeam is designed to process high time resolution (HTR) data recorded by the Voltage Capture System (VCS) subsystem of the Murchison Widefield Array (MWA) telescope.
It consists of

 1. A C library that can be used to create custom applications for processing VCS data,
 2. A suite of programs (see [Applications](@ref installationapplications) below).

## Dependencies

 - CUDA (required)
 - MPI
 - [PAL](https://github.com/Starlink/pal)
 - [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
 - [psrfits\_utils](https://github.com/demorest/psrfits_utils)
 - [mwa\_hyperbeam](https://github.com/mwatelescope/mwa_hyperbeam) (>=0.4.0)
 - [mwalib](https://github.com/MWATelescope/mwalib) (required)
 - [vdifio](https://github.com/demorest/vdifio)
 - [xGPU](https://github.com/GPU-correlators/xGPU)

## Compiling

cmake will attempt to locate the dependencies automatically. Assuming it succeeds, a typical (and minimal) cmake command will be

```bash
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=[target installation directory] \
      -DCMAKE_CUDA_ARCHITECTURES=[cuda_compute] \
      -DHYPERBEAM_HDF5=[path to mwa_full_embedded_element_pattern.h5] \
      ..
```

The HYPERBEAM\_HDF5 file can be supplied upon request, if it is not provided as part of the mwa\_hyperbeam library.

If some of the dependencies are in non-standard locations, cmake can be helped by setting the following cmake variables (using the -DVARIABLE=value syntax):
```bash
PAL_ROOT_DIR
CFITSIO_ROOT_DIR
PSRFITS_UTILS_ROOT_DIR
HYPERBEAM_ROOT
MWALIB_ROOT
VDIFIO_ROOT
```

Other optional variables can be set:
```bash
RUNTIME_DIR  -- Where to install needed runtime files (e.g. pq_phase_correction.txt)
```

Compilation and installation:
```bash
make
make install
```

## Applications {#installationapplications}

The following applications are built along with the vcsbeam library. The table lists each application's dependencies ('Y' = required, 'C' = required only at compile time):

|                                   | CUDA | MPI | PAL | cfitsio | psrfits\_utils | mwa\_hyperbeam | mwalib | vdifio | xGPU |
|-----------------------------------|:----:|:---:|:---:|:-------:|:--------------:|:--------------:|:------:|:------:|:----:|
| `fine_pfb_offline`                |   Y  |     |     |         |                |                |    Y   |        |      |
| `make_mwa_incoh_beam`             |   Y  |  Y  |     |    Y    |        Y       |                |    Y   |        |      |
| `make_mwa_tied_array_beam`        |   Y  |  Y  |  Y  |    Y    |        Y       |        Y       |    Y   |    Y   |      |
| `mwa_track_primary_beam_response` |   C  |     |  Y  |         |                |        Y       |    Y   |        |      |
| `mwa_mwa_tied_array_beam_psf`     |   C  |     |  Y  |         |                |        Y       |    Y   |        |      |
| `offline_correlator`              |   Y  |     |     |    Y    |                |                |    C   |        |   Y  |

Only those applications use dependencies are all present will be compiled and installed.
