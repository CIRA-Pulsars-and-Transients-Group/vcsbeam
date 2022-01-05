# VCSBeam Documentation

This documentation describes VCSBeam, a software package designed for processing high time resolution voltage data from the Murchison Widefield Array telescope.

## Source code

The source code for this software can be found [here](https://github.com/CIRA-Pulsars-and-Transients-Group/vcsbeam).

## Contents

 - Mathematical description
   + [Definitions](@ref definitions) -- Notation, coordinate systems
   + [Calibration](@ref calibration) -- RTS, Hyperdrive, Hyperbeam
   + [Beamforming](@ref beamforming) -- Applying calibration, phasing up the array, summing antennas
 - User's Guide
   + [Overview](@ref usersguideoverview)
   + Workflows
     * [Beamforming](@ref usersguidebeamforming)
     * [Preparing a calibration solution](@ref usersguidecalibration).
     * [Examples](@ref usersguideexamples)
   + Applications
     * [Fine PFB Offline](@ref applicationsfinepfboffline) -- Convert coarse channels to fine channels
     * [Make MWA Incoh Beam](@ref applicationsmakemwaincohbeam) -- Form an incoherent beam
     * [Make MWA Tied Array Beam](@ref applicationsmakemwatiedarraybeam) -- Form a tied-array beam
     * [MWA Plot Calibration](@ref applicationsmwaplotcalibration) -- Generate a plot of the calibration solution
     * [MWA Tied Array Beam PSF](@ref applicationsmwatiedarraybeampsf) -- Generate a plot of the tied array beam point spread function
     * [MWA Track Primary Beam Response](@ref applicationsmwatrackprimarybeamresponse) -- Track the sensitivity of a tied-array beam for a given RA/Dec as it passes through the MWA's primary beam
 - Appendices
   + [File formats](@ref fileformats) - RTS, Offringa
