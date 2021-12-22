# Calibration

[TOC]

Calibration is the process by which the temporal variations in the instrumental and/or ionospheric response to an incoming signal are characterised, in order to account for them when processing observational data.
A *calibration solution* is a mathematical operation that can be applied either to measured visibilities or raw voltages to correct for these variations, and recover the visibilities (or raw voltages) that would have been measured under ideal conditions.

For the MWA, calibration solutions are modelled as a tile-dependent, linear transformation of two orthogonal polarisations.
Thus, the full calibration solution for each tile is represented by a \f$2\times2\f$ Jones matrix, \f${\bf J}\f$, which is multiplied to the incident electric field, \f${\bf e}\f$, to predict the measured voltages, \f${\bf v}\f$.
\f[
    {\bf v} = {\bf J}{\bf e}.
\f]
The form of \f${\bf J}\f$ depends on the coordinate bases used for \f${\bf v}\f$ and \f${\bf e}\f$.
For example,
\f[
    \begin{bmatrix} v_q \\ v_p \end{bmatrix}
        = \begin{bmatrix}
            J_{qx} & J_{qy} \\
            J_{px} & J_{py}
        \end{bmatrix}
        \begin{bmatrix} e_x \\ e_y \end{bmatrix},
\f]
In practise, of course, the inverse operations are used to reconstruct the electric field from the measured voltages:
\f[
    {\bf e} = {\bf J}^{-1}{\bf v}.
\f]

Following [Sokolowski et al., 2017](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/calibration-and-stokes-imaging-with-full-embedded-element-primary-beam-model-for-the-murchison-widefield-array/FBA84B9EB94000BD6258A8F75840C476#), \f${\bf J}\f$ is calculated as the matrix product of two different matrices,
\begin{equation}
    {\bf J} = {\bf D}{\bf B},
\end{equation}
where
 - \f${\bf D}\f$ represents the (direction independent) instrumental gains, and
 - \f${\bf B}\f$ represents the (direction dependent) tile beam response.
\f${\bf D}\f$ is obtained either from the Real Time System, or Hyperdrive.
\f${\bf B}\f$ is obtained from Hyperbeam.

## Real Time System (RTS) -- [deprecated] {#rts}

The RTS is one of the pieces of software that can be used to generate calibration solutions for a VCS observation.
The solutions are given in two sets of files, called "DI_JonesMatrices_nodeCCC.dat" (hereafter, DIJones) and ``BandpassCalibration_nodeCCC.dat'' (hereafter, Bandpass), where CCC represents a coarse channel index.
The DIJones files contain both the Jones matrices for each tile for the given coarse channel, \f${\bf J}_d\f$, as well as the beam response matrix that was used by the RTS during calibration, \f${\bf B}_\text{cal}\f$.
To get the calibration solutions for the fine channels, one must multiply \f${\bf J}_d\f$ on the right by the matrices in the Bandpass files, \f${\bf J}_b\f$.

The bases of these matrices are as follows:
\f{align*}{
    {\bf J}_d &= \begin{bmatrix} J_{px} & J_{py} \\ J_{qx} & J_{qy} \end{bmatrix}, &
    {\bf J}_b &= \begin{bmatrix} J_{px} & J_{py} \\ J_{qx} & J_{qy} \end{bmatrix}, &
    {\bf B}_\text{cal} &= \begin{bmatrix} B_{px} & B_{py} \\ B_{qx} & B_{qy} \end{bmatrix}.
\f}
If you want to apply these calibration solutions to the same observation, then you can set
\f[
    {\bf J} = {\bf J}_d
\f]
for coarse channel approximations, or, for the fine channel solutions,
\f[
    {\bf J} = {\bf J}_d {\bf J}_b
\f]
If the calibration solution is to be applied to a different observation (in a different part of the sky), or even to a separate different direction in the same observation, then the instrumental gains must be reconstructed by factoring out the beam response:
\f[
    {\bf D} = {\bf J}{\bf B}_\text{cal}^{-1}.
\f]
The resulting instrumental gains matrix is in the basis
\f[
    {\bf D} = \begin{bmatrix} D_{pp} & D_{pq} \\ D_{qp} & D_{qq} \end{bmatrix}.
\f]

## Hyperdrive {#hyperdrive}

Hyperdrive outputs the instrumental gains matrices (per tile per fine channel) in the "Offringa" format, described in [Offringa format](#offringa).
The matrices are in the basis:
\f[
    {\bf D} = \begin{bmatrix}
        D_{qq} & D_{qp} \\
        D_{pq} & D_{pp}
    \end{bmatrix}.
\f]

## Hyperbeam {#hyperbeam}

The primary beam is obtained using [Hyperbeam](https://github.com/MWATelescope/mwa_hyperbeam), which is based on the FEE beam model described in [Sokolowski et al., 2017](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/calibration-and-stokes-imaging-with-full-embedded-element-primary-beam-model-for-the-murchison-widefield-array/FBA84B9EB94000BD6258A8F75840C476#).
This gives the \f${\bf B}\f$ matrix in the basis
\f[
    {\bf B}_\text{hb} = \begin{bmatrix} B_{q\theta} & B_{q\phi} \\ B_{p\theta} & B_{p\phi} \end{bmatrix}.
\f]
To obtain the beam matrix that converts from the \f$(x,y)\f$ basis to the \f$(q,p)\f$ basis, one must multiply \f${\bf B}_\text{hb}\f$ on the right by the parallactic angle correction matrix:
\f[
    {\bf B} = {\bf B}_\text{hb}{\bf P}_\text{pa}
            = \begin{bmatrix} B_{q\theta} & B_{q\phi} \\ B_{p\theta} & B_{p\phi} \end{bmatrix}
              \begin{bmatrix} P_{\theta x} & B_{\theta y} \\ B_{\phi x} & B_{\phi y} \end{bmatrix}
            = \begin{bmatrix} B_{qx} & B_{qy} \\ B_{px} & B_{py} \end{bmatrix}
\f]

