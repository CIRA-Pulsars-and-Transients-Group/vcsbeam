# User's Guide -- Beamforming {#usersguidebeamforming}

[TOC]

## Overview

To form a tied-array beam with a VCS observation, you will need both

 1. the data of the target observation, and
 2. an MWA calibration solution valid for the time period around which the target observation was made.

The following diagram gives an overview of the beamforming pipeline:

\dotfile beamforming.dot

Green nodes represent data files that are stored on disk, until deleted by the user.
Grey boxes represent applications, some of which are provided by VCSBeam, and some of which are external applications.
The links in the diagram take to you each application's documentation.

## Downloading the data

The new ASVO system for downloading is almost, but not quite ready for general use.
The documentation can be found [here](https://wiki.mwatelescope.org/display/MP/Data+Access).
For downloading Legacy (pre- September 2021) data, instructions can be found [here](https://wiki.mwatelescope.org/display/MP/Documentation#Documentation-Downloadingdatadownloading).
At some point, the deprecated Legacy method will no longer work, but when this will happen is currently unclear.

Note that the Legacy downloading method includes a built-in call to `recombine`, so that "downloading the data" will always result in a set of "Recombined voltages (.dat)" files.
However, there is currently no automated way of running `recombine` on Legacy data downloaded with the new ASVO system, and must be done manually.

Summary table for downloading instructions:
| Legacy | MWAX |
| ------ | ---- |
| [Documentation](https://wiki.mwatelescope.org/display/MP/Documentation#Documentation-Downloadingdatadownloading) | [Documentation](https://wiki.mwatelescope.org/display/MP/Data+Access) |

## Obtaining a calibration solution {#usersguidecalibration}

Which calibration software to use (RTS vs Hyperdrive) depends on whether the data is Legacy or MWAX, and whether it is VCS or already-correlated data (e.g. a dedicated calibration observation).
The following table summarises the possibilities:

[RTS]: https://wiki.mwatelescope.org/display/MP/Documentation#Documentation-CalibratingwiththeRealTimeSystem(RTS)
[Hyperdrive]: https://wiki.mwatelescope.org/pages/viewpage.action?pageId=52068764

|        | VCS                      | Correlated visibilities              |
| ------ | ------------------------ | ------------------------------------ |
| Legacy | [RTS][RTS]               | [RTS][RTS], [Hyperdrive][Hyperdrive] |
| MWAX   | [Hyperdrive][Hyperdrive] | [Hyperdrive][Hyperdrive]             |

The links in the table will take you to the corresponding documentation.
Apart from [the MWA Telescope Wiki][Hyperdrive] (same link as given in the table), Hyperdrive also has some documentation on [its Github main page](https://github.com/MWATelescope/mwa_hyperdrive), and [its Github Wiki page](https://github.com/MWATelescope/mwa_hyperdrive/wiki).

## Beamforming

### Set up
Once the data for the target observation have been downloaded, and a calibration solution obtained, the data may be beamformed.
This is achieved by using the application `make_mwa_tied_array_beam`.

`make_mwa_tied_array_beam` is highly configurable, and the full set of options can be seen by running
```
make_mwa_tied_array_beam -h
```
A minimal call requires
 1. The metafits file for the target observation,
 2. The metafits file for the calibration observation,
 3. The path to the calibration solution,
 4. A "pointings" file containing a list of RA/Decs.

#### Obtaining metafits files

The metafits files can be obtained via an [MWA Web Service](https://wiki.mwatelescope.org/display/MP/Web+Services), e.g.
```
wget -O <obsid>_metafits.fits http://ws.mwatelescope.org/metadata/fits?obs_id=<obsid>
```
where `<obsid>` is a placeholder for the MWA Observation ID.
For beamforming purposes, the metafits files can be given any arbitrary names.

#### The calibration path

The path to the calibration solution will be either the calibration file solution itself if it is an [Offringa-style format binary file](@ref offringafileformat), or the directory containing the calibration solution files if it is an [RTS solution](@ref rtsfileformat).

#### The pointings file

The "pointings" file is a text file containing one or more whitespace-separated RA/Dec pairs in the format
```
HH:MM:SS.S  DD:MM:SS.S
HH:MM:SS.S  DD:MM:SS.S
...
```
Eventually, pulsar names may also be used instead of RAs and Decs, but this is not yet supported.
To find the RA and Dec for one or more specific pulsars, use [the ATNF database](https://www.atnf.csiro.au/research/pulsar/psrcat/), or equivalently, the `psrcat` utility, e.g.:
```
psrcat -c "raj decj" -x B0031-07 J0437-4715 | awk '{print $1, $4}' > pointings.txt
```

### Choosing the output format

VCSBeam currently supports two output file formats:

 1. PSRFITS (full Stokes), 10 kHz frequency channels
 2. VDIF (veamformed voltages), 1.28 MHz frequency channels

Whether you choose PSRFITS or VDIF depends on your particular science case.
Since PSRFITS records **detected** samples, they cannot subsequently be coherently beamformed.
Therefore, S/N may be degraded if dispersion smearing across a 10 kHz channel is larger than the intrinsic pulse width.
For example, at 150 MHz, the dispersion smearing of PSR J2241-5236 across a 10 kHz channel is 0.28 ms, whereas the width of its pulse window is W50 = 0.07 ms, so its profile will be significantly smeared (see [Kaur et al., 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...882..133K/abstract)).
In contrast, VDIF files are coherently de-dispersed when processed with [DSPSR](http://dspsr.sourceforge.net/) and [PSRCHIVE](http://psrchive.sourceforge.net/).
The advantage to using PSRFITS is that PSRFITS is a supported format of [PRESTO](https://github.com/scottransom/presto), which is designed with searching for accelerated binary pulsars in mind.
The SMART survey uses PRESTO as part of its search pipeline.

VCSBeam can output either PSRFITS or VDIF formats, or both, by including the `-p` (PSRFITS) and `-v` (VDIF) options.
If neither `-p` nor `-v` is explicitly given, the default behaviour is to output only PSRFITS.

#### Channelisation options

Internally, VCSBeam always performs the beamforming operation on fine (10 kHz) channels.
Therefore, conversion to fine channels is necessary if the input data are MWAX coarse channels.
Conversion back to coarse channels is performed if VDIF output is requested.

\dotfile channlisation.dot

