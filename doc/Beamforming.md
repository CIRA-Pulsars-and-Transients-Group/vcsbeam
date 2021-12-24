# Beamforming {#beamforming}

[TOC]

Beamforming is expressed mathematically by the expression
\f[
    {\bf e}_f = \frac{1}{N_a} \sum_a e^{i\varphi_{a,f}} {\bf J}_{a,f}^{-1} {\bf v}_{a,f},
\f]
where

 - \f$a\f$ represents tiles (i.e. "antennas")
 - \f$f\f$ represents frequency channels
 - \f$\varphi\f$ is the delay phase,
 - \f${\bf v}\f$ is the Jones vector describing the instrumental voltages
 - \f${\bf J}\f$ is the Jones matrix describing the instrumental response to an incident electric field
 - \f${\bf e}\f$ is the Jones vector describing the (recovered) incident electric field
 - \f$N_a\f$ is the number of antennas

(The time dependence of these quantities is not shown here.)

Conceptually, we describe beamforming as consisting of the following distinct steps:

  1. Generating and applying the calibration Jones matrices to the voltages measured at each tile, (i.e. multiplying \f${\bf J}_{a,f}^{-1}\f$ to \f${\bf v}_{a,f}\f$),
  2. Shifting the signals from each tile in time to account for the signal delay from the look-direction to each tile, due to the geometric layout of the MWA tiles ("phasing up"), (i.e. multiplying \f$e^{i\varphi_{a,f}})\f$, and
  3. Averaging over the tiles.

For some operating modes, VCSBeam also performs *detection*, which means forming Stokes parameters, also described below.

These steps are described more fully in [Ord et al., 2019](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/abs/mwa-tiedarray-processing-i-calibration-and-beamformation/E9A7A9981AE9A935C9E08500CA6A1C1E), but are discussed further in the following subsections.

## Generating and applying the calibration solutions

The Jones matrix consists of a *direction independent* component, \f${\bf D}_a\f$, and a *direction dependent* component, \f${\bf B}_{a,f}\f$, such that
\f[
    {\bf J}_{a,f} = {\bf D}_a {\bf B}_{a,f}
\f]
\f${\bf D}_a\f$, which only depends on the antenna, is the *instrumental calibration solution*.
It is obtained using separate, dedicated calibration software, and VCSBeam currently supports solutions in two formats:

  1. The [RTS format](@ref rtsfileformat), and
  2. The [Offringa file format](@ref offringafileformat).

\f${\bf B}_{a,f}\f$, which formally depends on both antenna *and* frequency, is the *beam model*.
If all antennas were identical, the beam model would only depend on frequency, and we could write simply \f${\bf B}_f\f$.
In reality, antennas can differ because individual dipole elements can fail at different times.
These failures, when detected, are recorded in an observation's metadata, and are used by VCSBeam to obtain beam models for every configuration of live/dead dipoles that are present in a given observation.

The beam models themselves are calculated using [Hyperbeam](https://github.com/MWATelescope/mwa_hyperdrive), which implements the FEE beam model described in [Sokolowski et al. (2017)](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/calibration-and-stokes-imaging-with-full-embedded-element-primary-beam-model-for-the-murchison-widefield-array/FBA84B9EB94000BD6258A8F75840C476#).

The product obtained by multiplying \f${\bf J}^{-1}\f$ to the voltage data is
\f[
    \tilde{\bf e}_{a,f} = {\bf J}^{-1}_{a,f} {\bf v}_{a,f}.
\f]

## Phasing up

Phasing up refers to the process of accounting for the fact that each antenna will "see" an astrophysical signal arriving at a different time due to

  1. Its location relative to the source, and
  2. The different lengths of the cables that connect the antennas to the rest of the system.

Accounting for these time differences can be achieved either by a simple shift in the time domain, or, equivalently, the application of a phase ramp in the frequency domain, via the shift theorem.
The latter is what's implemented in VCSBeam, as indicated by the \f$e^{i\varphi}\f$ term in the expression above.
Thus, the result of phasing up is
\f[
    {\bf e}_{a,f} = e^{i\varphi_f} \tilde{\bf e}_{a,f}.
\f]

## Averaging the voltages

The final step is simply summing the voltages over all antennas.
As long as the correct delays have been applied, summing the voltages will coherently combine any astrophysical signal arriving from the specified look-direction.
The result is
\f[
    {\bf e}_f = \frac{1}{N_a}\sum_a \tilde{\bf e}_{a,f}.
\f]

## Detection (forming Stokes parameters)

VCSBeam converts the summed voltages into Stokes parameters if the PSRFITS output format is requested.
This is achieved by forming the *coherency matrix*:
\f[
    {\bf e}{\bf e}^\dagger
      = \begin{bmatrix}
            e_x e_x^\ast & e_x e_y^\ast \\
            e_y e_x^\ast & e_y e_y^\ast
        \end{bmatrix}
      = \frac12 \begin{bmatrix}
            I + Q & U + Vi \\
            U - Vi & I - Q
        \end{bmatrix}
\f]
where the frequency dependence of all terms is implicit.

However, as described in [Ord et al. (2019)](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/abs/mwa-tiedarray-processing-i-calibration-and-beamformation/E9A7A9981AE9A935C9E08500CA6A1C1E), VCSBeam also subtracts the autocorrelations, which for signals which are noise-dominated on a single-tile basis, improves the signal-to-noise ratio.
Thus, the actual detection operation that is implemented in VCSBeam is
\f[
    \frac12 \begin{bmatrix}
        I + Q & U + Vi \\
        U - Vi & I - Q
    \end{bmatrix}
    = {\bf e}{\bf e}^\dagger - \sum_a {\bf e}_a {\bf e}_a^\dagger.
\f]
