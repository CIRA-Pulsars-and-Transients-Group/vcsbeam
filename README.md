# VCSBeam

VCSBeam is designed to process high time resolution (HTR) data recorded by the Voltage Capture System (VCS) subsystem of the Murchison Widefield Array (MWA) telescope.
It consists of

 1. A C library that can be used to create custom applications for processing VCS data,
 2. A suite of programs (see [Applications](@ref installationapplications) below).

## Installation

See [the Installation Guide](https://cira-pulsars-and-transients-group.github.io/vcsbeam/installation.html).

## Alternatives

A Docker image containing (a possibly deprecated version of) similar beamforming code can be found [here](https://cloud.docker.com/u/cirapulsarsandtransients/repository/docker/cirapulsarsandtransients/vcstools).

## Documentation

 - [High-level description and examples of use on Garrawarla](https://wiki.mwatelescope.org/display/MP/Processing+high+time+resolution+data+with+VCSBeam) -- (deprecated: all documentation will eventually move to the Doxygen link below).
 - [Doxygen-generated reference documentation](https://cira-pulsars-and-transients-group.github.io/vcsbeam/)

## Credit

You can reference this repository using: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3762792.svg)](https://doi.org/10.5281/zenodo.3762792)

If you use the MWA beamformer, please give credit by citing:
[Ord et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019PASA...36...30O/abstract)

If you used polarimetry, please give credit by citing: 
[Xue et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019PASA...36...25X/abstract)

If you used the inverse PFB, please give credit by citing:
[McSweeney et al. (2020)](http://dx.doi.org/10.1017/pasa.2020.24)
