# User's Guide -- Beamforming {#usersguidebeamforming}

[TOC]

## Overview

To form a tied-array beam with a VCS observation, you will need both

 1. the data of the target observation, and
 2. an MWA calibration solution valid for the time period around which the target observation was made.

This page describes just the beamforming part (see [Calibration](@ref usersguidecalibration) for how to prepare a calibration solution).
The following diagram gives an overview of the beamforming pipeline:

\dotfile beamforming.dot

Green nodes represent data files that are stored on disk, until deleted by the user.
Down arrows represent applications, and those in blue are provided by VCSBeam.

## Downloading the data

The new ASVO system for downloading is almost, but not quite ready for general use.
The documentation can be found [here](https://wiki.mwatelescope.org/display/MP/Data+Access).
For downloading Legacy (pre- September 2021) data, instructions can be found [here](https://wiki.mwatelescope.org/display/MP/Documentation).
At some point, the deprecated Legacy method will no longer work, but when this will happen is currently unclear.

Note that the Legacy downloading method includes a built-in call to `recombine`, so that "downloading the data" will always result in a set of "Recombined voltages (.dat)" files.
However, there is currently no automated way of running `recombine` on Legacy data downloaded with the new ASVO system, and must be done manually.

## Obtaining a calibration solution

Which calibration software to use (RTS vs Hyperdrive) depends on whether the data is Legacy or MWAX, and whether you are doing in-beam calibration or using a dedicated calibration observation.
The following table summarises the possibilities:

|        | In-beam (voltages) | Dedicated calibration (calibrated voltages) |
| ------ | ------------------ | ------------------------------------------- |
| Legacy | RTS                | RTS, Hyperdrive                             |
| MWAX   | Hyperdrive         | Hyperdrive                                  |

### The Real Time System (RTS)

Documentation for using the RTS can be found with [the VCSTools documentation](https://wiki.mwatelescope.org/display/MP/Documentation).

### Hyperdrive

Documentation for Hyperdrive can be found on [this MWA Telescope Wiki page](https://wiki.mwatelescope.org/pages/viewpage.action?pageId=52068764), [its Github main page](https://github.com/MWATelescope/mwa_hyperdrive), and [its Github Wiki page](https://github.com/MWATelescope/mwa_hyperdrive/wiki).

## Beamforming

The main 
