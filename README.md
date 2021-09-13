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
```

Other optional variables can be set:
```bash
RUNTIME_DIR  -- Where to install needed runtime files (e.g. pq_phase_correction.txt)
```

Upon successful completion of the cmake command, the following will build and install make\_beam:
```bash
make
make install
```

A Docker image containing (a possibly deprecated version of) make\_beam can be found [here](https://cloud.docker.com/u/cirapulsarsandtransients/repository/docker/cirapulsarsandtransients/vcstools)

Help
------
Documentation on how to process MWA VCS data, including the use of make\_beam, can be found [here](https://wiki.mwatelescope.org/display/MP/Documentation)

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
