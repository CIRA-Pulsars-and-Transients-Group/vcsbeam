VCSBeam
======

Installation
------

### Dependencies

 - CUDA
 - [PAL](https://github.com/Starlink/pal)
 - [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
 - [psrfits\_utils](https://github.com/demorest/psrfits_utils)
 - [mwa\_hyperbeam](https://github.com/mwatelescope/mwa_hyperbeam)
 - [mwalib](https://github.com/MWATelescope/mwalib)
 - [vdifio](https://github.com/demorest/vdifio)
 - [xGPU](https://github.com/GPU-correlators/xGPU) (optional)

### Compiling

cmake will attempt to locate the dependencies automatically. Assuming it succeeds, a typical (and minimal) cmake command will be

```bash
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=[target installation directory] \
      -DCMAKE_CUDA_ARCHITECTURES=[cuda_compute] \
      ..
```

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

The Hyperbeam dependency also requires the file `mwa_full_embedded_element_pattern.h5` (which it supplies), the full path to which must be passed to cmake with the option
```
-DHYPERBEAM_HDF5=/path/to/mwa_full_embedded_element_pattern.h5
```

Upon successful completion of the cmake command, the following will build and install the vcsbeam library, as well as several applications that make use of it (see below):
```bash
make
make install
```

Exectuable Programs (Apps)
------
Currently the following apps are built along with the vcsbeam library:
1. `fine_pfb_offline`
2. `make_mwa_incoh_beam`
3. `make_mwa_tied_array_beam`
4. `mwa_track_primary_beam_response`

Their dependencies are as follows (however, currently all but xGPU are "required" at compile time):

|                                   | PAL | cfitsio | psrfits\_utils | mwa\_hyperbeam | mwalib | vdifio | xGPU |
|-----------------------------------|:---:|:-------:|:--------------:|:--------------:|:------:|:------:|:----:|
| `fine_pfb_offline`                |     |         |                |                |    Y   |        |      |
| `make_mwa_incoh_beam`             |     |    Y    |        Y       |                |    Y   |        |      |
| `make_mwa_tied_array_beam`        |  Y  |    Y    |        Y       |        Y       |    Y   |    Y   |      |
| `mwa_track_primary_beam_response` |  Y  |         |                |        Y       |    Y   |        |      |
| `offline_correlator`              |     |         |                |                |        |        |   Y  |

Alternatives
------
A Docker image containing (a possibly deprecated version of) similar beamforming code can be found [here](https://cloud.docker.com/u/cirapulsarsandtransients/repository/docker/cirapulsarsandtransients/vcstools)

Help
------
(Still coming!)

Credit
------
You can reference this repository using: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3762792.svg)](https://doi.org/10.5281/zenodo.3762792)

If you use the MWA beamformer, please give credit by citing:
[Ord et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019PASA...36...30O/abstract)

If you used polarimetry, please give credit by citing: 
[Xue et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019PASA...36...25X/abstract)

If you used the inverse PFB, please give credit by citing:
[McSweeney et al. (2020)](http://dx.doi.org/10.1017/pasa.2020.24)
