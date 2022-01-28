# User's Guide -- Beamforming {#usersguidebeamforming}

[TOC]

## Overview

To form a tied-array beam with a VCS observation, you will need both

 1. the data of the target observation, and
 2. an MWA calibration solution valid for the time period around which the target observation was made.

The following diagram gives an overview of the beamforming pipeline:

\dotfile beamforming.dot

Blue nodes represent data files that are stored on disk, until deleted by the user.
Grey boxes represent applications, some of which are provided by VCSBeam, and some of which are external applications.
The links in the diagram take to you each application's documentation.

________________

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

________________

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

________________

## Beamforming

________________

### Setting up
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

________________

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

________________

### Input options

The default behaviour of [make_mwa_tied_array_beam](@ref applicationsmakemwatiedarraybeam) is to look for the input files in the current working directory, and to process all time steps, and as many coarse channels as MPI processes are used.
If the input data are in another directory, this directory can be passed to the `-d` option.

#### Choosing which timesteps to process

You can set the range of timesteps to be processed via the `-b` and `-T` options.

The argument of `-b` is the beginning time, which can be either an absolute time (GPS second) or a relative time.
To indicate a relative time, the first character of the argument must be either '`+`' or '`-`'.
If it is '`+`', then the number is considered an offset from the first "good" second, where "good" is defined in the metafits file as the first second after the "quack time" has elapsed.
If it is '`-`', then the number is considered an offset from the *end* of the observation, with `-1` indicating the last second (similar to Python-style indexing).
The default value is `+0`.

The argument of `-T` is the number of seconds of data to process.
If not supplied, it will process all available seconds from the specified beginning time onwards.

\todo There might be a bug whereby using all the defaults crashes because the default total number of seconds to be processed is more than the number of seconds after the "good" time starts.

#### Choosing which coarse channels to process

You can set the range of coarse channels to be processed via the `-f` option and by setting the number of MPI processes, with one coarse channel processed per process.

The argument of `-f` is the lowest coarse channel to be processed.
It can either be an MWA receiver channel number (equal to the coarse channel centre frequency divided by 1.28 MHz), or a relative receiver channel number.
To indicate a relative channel number, the first character of the argument must be either '`+`' or '`-`'.
If it is '`+`', then the number is considered an offset from the lowest coarse channel.
If it is '`-`', then the number is considered an offset from the highest coarse channel plus one, with `-1` therefore indicating the highest channel (similar to Python-style indexing).
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

________________

### Calibration options

[make_mwa_tied_array_beam](@ref applicationsmakemwatiedarraybeam) will expect an RTS style solution unless the `-O` option is explicitly given (this default behaviour may change in the future as Hyperdrive eventually supercedes the RTS as the primary calibration tool used for VCS data).

The calibration solutions provided by the RTS or Hyperdrive can be further manipulated in several ways before they are applied to the input voltages.
Many of these manipulates are still experimental, and mostly affect the fidelity of the polarisation response, whose verification is still a work in progress.

#### Flagging extra tiles

Extra tiles can be flagged by passing a text file containing TileNames to be flagged to the `-F` option.
(To decide *which* tiles need flagging, and how to prepare this file, see [Preparing a calibration solution](@ref usersguidecalibration).)

#### Including the Bandpass solutions

The `-B` option signals that the fine channel calibration solutions should be used.
This only applies to the RTS solutions, where the solutions are separated out into "coarse channel" (DIJones) and "fine channel" (Bandpass) solutions (see [Calibration](@ref calibration) for more details).

The advantage of using the Bandpass solutions is that the solutions may be more accurate (i.e. more accurately reflect the true instrumental response) for individual fine channels, whereas the DIJones solutions, by themselves, only include a zeroth order (i.e. constant) approximation to the solution for the whole coarse channel.
If the phases do not not change with respect to frequency too rapidly, then the DIJones solutions are probably a good enough approximation and do not degrade the S/N too much.
If the phases *do* change rapidly for a given antenna, then that antenna's contribution near the coarse channel edges will be degraded to some extent.

Ironically, the edge channels are usually flagged (by default) when producing the calibration solution in the first place.
Thus, if the Bandpass solutions are used, the edge channels "remain" flagged during beamforming, so any signal present there will be lost anyway.

#### Keeping the cross terms?

By default, the off-diagonal ("cross") terms of the Jones matrices (PQ and QP) are set to zero, based on the premise that these terms are dominated by noise.
This claim is tantamount to saying that there is negligible instrumental leakage between the two polarisations (P and Q).
Anecdotally, the MWA imaging group claim that in their experience, the solutions produce better images when the cross terms are zeroed.

The `-X` option can be used to retain the cross terms.

#### Dividing a reference antenna

The calibration solution is a Jones matrix whose absolute phase carries no physical signifance.
Said another way, calibration solutions are unique only up to a complex unit scalar factor, so that if \f$\{{\bf J}_1, {\bf J}_2, \dots, {\bf J}_{N_a}\}\f$ is a set of calibration solutions for antennas \f$1, 2, \dots, N_a\f$, then \f$\{e^{i\theta}{\bf J}_1, e^{i\theta}{\bf J}_2, \dots, e^{i\theta}{\bf J}_{N_a}\}\f$ is an equivalent set of solutions, for any arbitrary \f$\theta\f$.
Thus, the solutions produced by the RTS and/or Hyperdrive may "look" worse than they really are because the solutions for different frequency channels may have different arbitrary phases.
Thus, in order to compare the solutions across frequency, it is useful to divide the elements of all Jones matrix by a unit complex number which has the same phase as an arbitrarily chosen reference antenna.
This will make the reference antenna itself have a calibration solution consisting only of real-valued elements, but the phases of the other antennas can now be more easily evaluated for "goodness".

\todo Finish the calibration section (and probably move it all to the Calibration workflow page).

________________

### Memory management

With such large data sizes, it is possible that there is not enough memory on the GPU to carry all the necessary data arrays to process even one second of one coarse channel!
(This is almost certainly true for desktop computers, but probably not true for supercomputers.)

Currently, the amount of memory needed is not calculated (but will be implemented in future).
If the required memory is larger than the available memory, this problem can be mitigated by using the `-n` option, which tells [make_mwa_tied_array_beam](@ref applicationsmakemwatiedarraybeam) to process only 1/`nchunks` seconds of data at a time, where `nchunks` is the argument to `-n`.
If `nchunks` does not divide evenly into the number of samples per second of fine channelised data (10,000), then it is automatically increased until it does.

________________

## Examples

### MWAX observation of PSR B0031-07

 1. Get the metafits files for the target observation and the calibration observation:
```
wget -O 1320499816.fits http://ws.mwatelescope.org/metadata/fits?obs_id=1320499816
wget -O 1320412440.fits http://ws.mwatelescope.org/metadata/fits?obs_id=1320412440
```

 2. Generate a pointing file for B0031-07:
```
psrcat -c "raj decj" -x B0031-07 | awk '{print $1, $4}' > pointings.txt
```
which produces the file `pointings.txt`:
```
00:34:08.8703 -07:21:53.409
```

 3. The calibration solution was obtained by the method given in the example on the [Calibration](@ref usersguidecalibration) page.
During that process, it was discovered that tile HexE2 needed to be flagged, so a file `flagged_tilenames.txt` was created with the following single entry:
```
HexE2
```
Otherwise, the solutions already looked good enough, so that no further amendments were deemed necessary.

 4. After deciding how much of the data to process, and what I wanted for the output format, the final [make_mwa_tied_array_beam](@ref applicationsmakemwatiedarraybeam) command is:
```
mpirun -np 24 make_mwa_tied_array_beam \      # Use 24 MPI processes to beamform 24 coarse channels and splice them together in the output PSRFITS
    -m 1320499816.fits \                      # The target observation metafits file
    -c 1320412440.fits \                      # The calibration observation metafits file
    -C 1320412440_hyperdrive_solutions.bin \  # The calibration solution
    -O \                                      # Signal that the above calibration solution is an Offringa-style solution
    -d /astro/mwavcs/asvo/252057 \            # The location of the MWAX voltage (.sub) files
    -P pointings.txt \                        # The pointings file
    -b 1320499824 \                           # Start 8 seconds in
    -T 592 \                                  # Process 592 seconds
    -f 109 \                                  # Start at channel 109
    -F flagged_tilenames.txt \                # The flagged channels file
    -p \                                      # Output PSRFITS (10 kHz, full Stokes)
    -R NONE \                                 # Do not divide through any reference tile
    -U 0,0 \                                  # Do not apply any extra phase ramp
    -X                                        # Keep the cross terms
```
