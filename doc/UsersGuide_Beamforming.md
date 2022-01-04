# User's Guide -- Beamforming {#usersguidebeamforming}

[TOC]

## Overview

To form a tied-array beam with a VCS observation, you will need both

 1. the data of the target observation, and
 2. an MWA calibration solution valid for the time period around which the target observation was made.

This page describes just the beamforming part (see [Calibration](@ref usersguidecalibration) for how to prepare a calibration solution).
An overview of the entire beamforming pipeline is represented in the following diagram:

\dotfile beamforming.dot

## Downloading the data

The new system for downloading is almost, but not quite ready for general use.
The documentation can be found [here](https://wiki.mwatelescope.org/display/MP/Data+Access).
For downloading Legacy (pre- September 2021) data, instructions can be found [here](https://wiki.mwatelescope.org/display/MP/Documentation).
At some point, the deprecated Legacy method will no longer work, but when this will happen is currently unclear.

## Beamforming

The main
