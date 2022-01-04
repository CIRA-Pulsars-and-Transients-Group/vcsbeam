# User's Guide -- Beamforming {#usersguidebeamforming}

[TOC]

## Overview

To form a tied-array beam with a VCS observation, you will need both

 1. the data of the target observation, and
 2. an MWA calibration solution valid for the time period around which the target observation was made.

The following diagram gives an overview of the beamforming pipeline:

\dotfile beamforming.dot

Green nodes represent data files that are stored on disk, until deleted by the user.
Down arrows represent applications, and those in blue are provided by VCSBeam.

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

The main 
