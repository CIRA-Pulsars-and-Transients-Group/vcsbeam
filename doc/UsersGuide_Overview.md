# User's Guide -- Overview {#usersguideoverview}

[TOC]

@todo Add link to program names once those pages are written.

# VCSBeam Overview

VCSBeam is a software library and a suite of applications designed to process high time resolution (HTR) data recorded by the Voltage Capture System (VCS) subsystem of the Murchison Widefield Array (MWA) telescope.
The software repository for VCSBeam is located at [this GitHub page](https://github.com/CIRA-Pulsars-and-Transients-Group/vcsbeam).

The primary goal of VCSBeam is to combine the signals measured at each of the MWA's antenna elements ("tiles") in order to maximise the sensitivity of the entire array towards a desired look-direction.
(The beam pattern of the combined array is called a "tied-array beam", and the processing of forming it, "beamforming".)
VCSBeam provides the executable `make_mwa_tied_array_beam` to perform this beamforming operation.

VCSBeam includes a variety of tools to process HTR data in many ways, including (but not limited to):
 - Forming tied-array beams
 - Forming incoherent beams
 - Applying a polyphase filterbank (PFB)
 - Visualising calibration solutions
 - Visualising the tied-array beam pattern
 - Modelling the tied-array beam response towards a particular target

Each of the above items corresponds to an "application" that comes shipped with VCSBeam, and clicking the links will take you to each application's documentation.
However, these processing tasks are not the only scientifically useful ways that VCS data can be processed, and the broader goal of VCSBeam is to provide a flexible framework that can be used to implement other novel algorithms and processing techniques.
As such, VCSBeam provides a C library which can be used to create custom processing applications.
This library is still a work-in-progress, and documentation for its use will be added in the future.
