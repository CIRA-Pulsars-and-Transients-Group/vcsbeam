# User's Guide -- Overview {#usersguideoverview}

[TOC]

@todo Add link to program names once those pages are written.

# VCSBeam Overview

VCSBeam is a software library and a suite of applications designed to process high time resolution (HTR) data recorded by the Voltage Capture System (VCS) subsystem of the Murchison Widefield Array (MWA) telescope.
The software repository for VCSBeam is located at [this GitHub page](https://github.com/CIRA-Pulsars-and-Transients-Group/vcsbeam).

The primary goal of VCSBeam is to combine the signals measured at each of the MWA's antenna elements ("tiles") in order to maximise the sensitivity of the entire array towards a desired look-direction.
(The beam pattern of the combined array is called a "tied-array beam", and the processing of forming it, "beamforming".)
Processing an observation in this way is described on [this page](@ref usersguidebeamforming).

VCSBeam can not only beamform, but also process VCS data in other ways.
These include:

 - [Forming tied-array beams](@ref applicationsmakemwatiedarraybeam)
 - [Forming incoherent beams](@ref applicationsmakemwaincohbeam)
 - [Applying a polyphase filterbank (PFB) for fine channelisation](@ref applicationsfinepfboffline)
 - [Visualising calibration solutions](@ref applicationsmwaplotcalibration)
 - [Visualising the tied-array beam pattern](@ref applicationsmwatiedarraybeampsf)
 - [Modelling the tied-array beam response towards a particular target](@ref applicationsmwatrackprimarybeamresponse)

Each of the above items corresponds to an "application" that comes shipped with VCSBeam, and clicking the links will take you to each application's documentation.
However, these processing tasks are not the only scientifically useful ways that VCS data can be processed, and the broader goal of VCSBeam is to provide a flexible framework that can be used to implement other novel algorithms and processing techniques.
As such, VCSBeam provides a C library which can be used to create custom processing applications.
This library is still a work-in-progress, and documentation for its use will be added in the future.
