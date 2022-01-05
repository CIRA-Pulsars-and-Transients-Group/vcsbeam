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

## Obtaining a calibration solution

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

The [RTS][RTS] link describes a workflow for preparing a calibration solution using the RTS.
The equivalent workflow for Hyperdrive is found at the page [Preparing a calibration solution](@ref usersguidecalibration).
However, it should be noted that the visualisation tools used for Hyperdrive can also be used for RTS solutions.

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

### Output options

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
If neither `-p` nor `-v` is explicitly given, the default behaviour is to output the same channelisation as the input (i.e. MWAX&rarr;VDIF, Legacy&rarr;PSRFITS).

By default, PSRFITS files are written out with 200 seconds per file.
This behaviour can be altered with the `-t` option.

#### Channelisation options

Currently, only coarse (1.28 MHz) channels and fine (10 kHz) channels are supported, but this will eventually be generalised to allow arbitrary channelisation.

Internally, VCSBeam always performs the beamforming operation on fine (10 kHz) channels.
Therefore, conversion to fine channels is necessary if the input data are MWAX coarse channels.
Conversion back to coarse channels is performed if VDIF output is requested.
This is illustrated in the following diagram (as well as in table form in the docstring for the function vmSetOutputChannelisation()):

\dotfile channelisation.dot

VCSBeam allows the filters used for the analysis and synthesis PFBs to be chosen by the user, via the options `-A` and `-S`, respectively.
The available filters are supplied in [the pfb_filter folder](https://github.com/CIRA-Pulsars-and-Transients-Group/vcsbeam/tree/main/pfb_filter), which are copied to `RUNTIME_DIR` during installation.
If no filters are explicitly given, the default filters are used: `FINEPFB` for the analysis filter, `LSQ12` for the synthesis filter.
`FINEPFB` is the same filter that was implemented on Legacy MWA FPGAs.

Custom filters can also be used by placing the filter coefficients in a file with the name `<FILTER_NAME>.dat` and placing the file in `RUNTIME_DIR`.
Only analysis and synthesis filters with 12 taps (and therefore 12x128 = 1536 coefficients) have been tested, and behaviour for filters with different lengths is currently undefined.
Additionally, [make_mwa_tied_array_beam](@ref applicationsmakemwatiedarraybeam) currently only supports critically-sampled PFBs, although this will also be generalised in the future.

**WARNING**: the analysis filter coefficients **MUST** be integers.
A scaling factor is applied in order to ensure that the output fits the VDIF dynamic range.
Therefore, if a custom analysis filter is used, it is highly recommended to scale the coefficients so that (1) they are integers, and (2) they are normalised in the same way as `FINEPFB` (i.e. they sum to 15106424).

##### Fine channel sample formats (SMART survey option)

[McSweeney2020]: https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/mwa-tiedarray-processing-iii-microsecond-time-resolution-via-a-polyphase-synthesis-filter/C76837FE84A1DCFAA696C9AC10C44F2D

The Legacy MWA pdemoted the output of the analysis PFB to (4+4)bit complex integers. While this saves space, it marginally degrades the precision of the beamformer output.
As the fine channelisation of MWAX data is now done entirely in (GPU) memory, the demotion is no longer necessary; however, to ensure that the SMART survey data (which spans both Legacy and MWAX phases) is processed in a homogeneous way, an option to replicate the same demotion step prior to beamforming, the `-s` option is provided.
It should be noted that this demotion step implements an "asymmetric rounding" scheme, which is described in the Appendix of [McSweeney et al. (2020)][McSweeney2020].
Therefore, the use of this option is not recommended, outside of the express purpose of producing fine channels equivalent to MWA Phase 1 & 2.

### Input options

The default behaviour of [make_mwa_tied_array_beam](@ref applicationsmakemwatiedarraybeam) is to look for the input files in the current working directory, and to process all time steps, and as many coarse channels as MPI processes are used.
If the input data are in another directory, this directory can be passed to the `-d` option.

#### Choosing which timesteps to process

You can set the range of timesteps to be processed via the `-b` and `-T` options.

`-b` is the beginning time, and can be either an absolute time (GPS second) or a relative time.
To indicate a relative time, the first character of the argument must be either '`+`' or '`-`'.
If it is `'+'`, then the number is considered an offset from the first "good" second, where "good" is defined in the metafits file as the first second after the "quack time" has elapsed.
If it is `'-'`, then the number is considered an offset from the *end* of the observation, with `-1` indicating the last second (similar to Python-style indexing).
The default value for `-b` is `+0`.

`-T` is the number of seconds of data to process.
If not supplied, it will process all available seconds from the specified beginning time onwards.

\todo There might be a bug whereby using all the defaults crashes because the default total number of seconds to be processed is more than the number of seconds after the "good" time starts.

#### Choosing which coarse channels to process

You can set the range of coarse channels to be processed via the `-f` option and by setting the number of MPI processes, with one coarse channel processed per process.

`-f` is the lowest coarse channel to be processed.
It can either be an MWA receiver channel number (equal to the coarse channel centre frequency divided by 1.28 MHz), or a relative receiver channel number.
To indicate a relative channel number, the first character of the argument must be either '`+`' or '`-`'.
If it is `'+'`, then the number is considered an offset from the lowest coarse channel.
If it is `'-'`, then the number is considered an offset from the highest coarse channel plus one, with `-1` therefore indicating the highest channel (similar to Python-style indexing).
The default value is `+0`.

The number of channels processed is equal to the number of MPI processes chosen.
For example,
```
mpirun -np 6 make_mwa_tied_array_beam ...
```
will process the lowest 6 coarse channels in the observation,
```
mpirun -np 6 make_mwa_tied_array_beam -f -6 ...
```
will process the highest 6 channels, and
```
mpirun -np 6 make_mwa_tied_array_beam -f 115 ...
```
will process channels 115, 116, 117, 118, 119, and 120.

For VDIF output, each coarse channel is written to a separate file, so there is no difference between running, e.g.
```
mpirun -np 24 make_mwa_tied_array_beam -v ...
```
and
```
mpirun -np 1 make_mwa_tied_array_beam -v -f +0 ...
mpirun -np 1 make_mwa_tied_array_beam -v -f +1 ...
mpirun -np 1 make_mwa_tied_array_beam -v -f +2 ...
...
```
However, for PSRFITS output, the coarse channels are spliced together before being written out.
Therefore, running
```
mpirun -np 24 make_mwa_tied_array_beam -p ...
```
will result in a single PSRFITS file with all 24 coarse channels included, whereas running
```
mpirun -np 1 make_mwa_tied_array_beam -p -f +0 ...
mpirun -np 1 make_mwa_tied_array_beam -p -f +1 ...
mpirun -np 1 make_mwa_tied_array_beam -p -f +2 ...
...
```
will result in 24 separate PSRFITS files, one per coarse channel.
The advantage to running a single multi-process MPI job is that splicing is done automatically.
The disadvantage is that this can potentially lead to slower wall time completion of the beamforming job, since the splicing occurs after each second of data, requiring that the MPI processes are synchronised after processing each second of data.
However, this is unlikely to affect the wall time significantly (although this remains to be benchmarked) since the splicing is performed synchronously with the reading in of the subsequent second's worth of data, which is believed to be the current bottleneck for throughput.

### Calibration options

[make_mwa_tied_array_beam](@ref applicationsmakemwatiedarraybeam) will expect an RTS style solution unless the `-O` option is explicitly given (this default behaviour may change in the future as Hyperdrive eventually supercedes the RTS as the primary calibration tool used for VCS data).

The calibration solutions provided by the RTS or Hyperdrive can be further manipulated in several ways before they are applied to the input voltages.
Many of these manipulates are still experimental, and mostly affect the fidelity of the polarisation response, whose verification is still a work in progress.



#### Flagging extra tiles

Extra tiles can be flagged by passing a text file containing TileNames to be flagged to the `-F` option.
(To decide *which* tiles need flagging, and how to prepare this file, see [Preparing a calibration solution](@ref usersguidecalibration).)

#### Including the Bandpass solutions

The `-B` option signals that the fine channel calibration solutions should be used.
This only applies to the RTS solutions, where the solutions are separated out into "coarse channel" (DIJones) and "fine channel" (Bandpass) solutions (see.
