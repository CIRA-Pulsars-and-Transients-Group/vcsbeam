# Calibration {#calibration}

[TOC]

Calibration is the process by which the temporal variations in the instrumental
and/or ionospheric response to an incoming signal are characterised, in order to
account for them when processing observational data.
A *calibration solution* is a mathematical operation that can be applied either
to measured visibilities or raw voltages to correct for these variations, and
recover the visibilities (or raw voltages) that would have been measured under
ideal conditions.

For the MWA, calibration solutions are modelled as a tile-dependent, linear
transformation of two orthogonal polarisations.
Thus, the full calibration solution for each tile is represented by a
\f$2\times2\f$ Jones matrix, \f${\bf J}\f$, which is multiplied to the incident
electric field, \f${\bf e}\f$, to predict the measured voltages, \f${\bf v}\f$.
\f[
    {\bf v} = {\bf J}{\bf e}.
\f]
The form of \f${\bf J}\f$ depends on the [coordinate
bases](@ref coordinatesystems) used for \f${\bf v}\f$ and \f${\bf e}\f$.
For example,
\f[
    \begin{bmatrix} v_q \\ v_p \end{bmatrix}
        = \begin{bmatrix}
            J_{qx} & J_{qy} \\
            J_{px} & J_{py}
        \end{bmatrix}
        \begin{bmatrix} e_x \\ e_y \end{bmatrix},
\f]
In practise, of course, the inverse operations are used to reconstruct the
electric field from the measured voltages:
\f[
    {\bf e} = {\bf J}^{-1}{\bf v}.
\f]

Following [Sokolowski et al.,
2017](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/calibration-and-stokes-imaging-with-full-embedded-element-primary-beam-model-for-the-murchison-widefield-array/FBA84B9EB94000BD6258A8F75840C476#),
\f${\bf J}\f$ is calculated as the matrix product of two different matrices,
\f[
    {\bf J} = {\bf D}{\bf B},
\f]
where

- \f${\bf D}\f$ represents the (direction independent) instrumental gains, and
- \f${\bf B}\f$ represents the (direction dependent) tile beam response.

\f${\bf D}\f$ is obtained from Hyperdrive.
\f${\bf B}\f$ is obtained from Hyperbeam.

## Obtaining the instrumental gains

### Hyperdrive {#hyperdrive}

Hyperdrive is used to generate the calibration solutions. VCSBeam reads its
native FITS solutions directly ([Hyperdrive FITS
format](@ref hyperdrivefileformat)); this is the default. The "Offringa" binary
format ([Offringa file format](@ref offringafileformat)), which hyperdrive can
convert its solutions to, is also supported via the `-O` option. Hyperdrive
writes one file per contiguous subband. Each coarse channel is matched to its
solution channels by frequency, so the files may cover anywhere from a single
coarse channel up to the full 24. Each solution is applied to all of the voltage
fine channels that fall within the same solution channel (nearest-neighbour
replication). Flagged tiles, flagged channels and NaN solutions are zeroed, i.e.
flagged in the beamformer. The matrices are in the basis:
\f[
    {\bf D} = \begin{bmatrix}
        D_{qq} & D_{qp} \\
        D_{pq} & D_{pp}
    \end{bmatrix}.
\f]

## Obtaining the primary beam model

### Hyperbeam {#hyperbeam}

The primary beam is obtained using
[Hyperbeam](https://github.com/MWATelescope/mwa_hyperbeam), which is based on
the FEE beam model described in [Sokolowski et al.,
2017](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/calibration-and-stokes-imaging-with-full-embedded-element-primary-beam-model-for-the-murchison-widefield-array/FBA84B9EB94000BD6258A8F75840C476#).
This gives the \f${\bf B}\f$ matrix in the basis
\f[
    {\bf B}_\text{hb} = \begin{bmatrix} B_{q\theta} & B_{q\phi} \\ B_{p\theta} & B_{p\phi} \end{bmatrix}.
\f]
To obtain the beam matrix that converts from the \f$(x,y)\f$ basis to the
\f$(q,p)\f$ basis, one must multiply \f${\bf B}_\text{hb}\f$ on the right by the
parallactic angle correction matrix:
\f[
    {\bf B} = {\bf B}_\text{hb}{\bf P}_\text{pa}
            = \begin{bmatrix} B_{q\theta} & B_{q\phi} \\ B_{p\theta} & B_{p\phi} \end{bmatrix}
              \begin{bmatrix} P_{\theta x} & P_{\theta y} \\ P_{\phi x} & P_{\phi y} \end{bmatrix}
            = \begin{bmatrix} B_{qx} & B_{qy} \\ B_{px} & B_{py} \end{bmatrix}
\f]

