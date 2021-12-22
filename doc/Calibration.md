# Calibration

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

## Real Time System (RTS) -- [deprecated]

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
\begin{equation}
    {\bf J} = {\bf J}_d
\end{equation}
for coarse channel approximations, or, for the fine channel solutions,
\begin{equation}
    {\bf J} = {\bf J}_d {\bf J}_b
\end{equation}
If the calibration solution is to be applied to a different observation (in a different part of the sky), or even to a separate different direction in the same observation, then the instrumental gains must be reconstructed by factoring out the beam response:
\begin{equation}
    {\bf D} = {\bf J}{\bf B}_\text{cal}^{-1}.
\end{equation}
The resulting instrumental gains matrix is in the basis
\begin{equation}
    {\bf D} = \begin{bmatrix} D_{pp} & D_{pq} \\ D_{qp} & D_{qq} \end{bmatrix}.
\end{equation}

\subsection{Hyperdrive}
\label{sec:hyperdrive}
\index{hyperdrive}

Hyperdrive outputs the instrumental gains matrices (per tile per fine channel) in the ``Offringa'' format, described in Appendix \ref{app:offringa}.
The matrices are in the basis:
\begin{equation}
    {\bf D} = \begin{bmatrix}
        D_{qq} & D_{qp} \\
        D_{pq} & D_{pp}
    \end{bmatrix}.
\end{equation}

\subsection{Hyperbeam}
\label{sec:hyperbeam}
\index{hyperbeam}

The primary beam is obtained using \href{https://github.com/MWATelescope/mwa_hyperbeam}{Hyperbeam}, which is based on the FEE beam model\index{primary beam!full embedded element (FEE)} \citep{Sokolowski2017}.
This gives the \f${\bf B}\f$ matrix in the basis
\begin{equation}
    {\bf B}_\text{hb} = \begin{bmatrix} B_{q\theta} & B_{q\phi} \\ B_{p\theta} & B_{p\phi} \end{bmatrix}.
\end{equation}
To obtain the beam matrix that converts from the \f$(x,y)\f$ basis to the \f$(q,p)\f$ basis, one must multiply \f${\bf B}_\text{hb}\f$ on the right by the parallactic angle correction matrix:
\begin{equation}
    {\bf B} = {\bf B}_\text{hb}\pamat
            = \begin{bmatrix} B_{q\theta} & B_{q\phi} \\ B_{p\theta} & B_{p\phi} \end{bmatrix}
              \begin{bmatrix} P_{\theta x} & B_{\theta y} \\ B_{\phi x} & B_{\phi y} \end{bmatrix}
            = \begin{bmatrix} B_{qx} & B_{qy} \\ B_{px} & B_{py} \end{bmatrix}
\end{equation}

\section{Beamforming}

Beamforming consists of three steps:
\begin{enumerate}
    \item Applying the calibration Jones matrices to the voltages measured at each tile,
    \item Shifting the signals from each tile in time to account for the signal delay from the look-direction to each tile, due to the geometric layout of the MWA tiles (``phasing up''),
    \item Summing the voltages from all the tiles together.
\end{enumerate}
These steps are described more fully in \citet{Ord2019}, but are discussed further in the following subsections.

\subsection{Applying the calibration solutions}

In \vcsbeam{}, the following bases are used:
\begin{equation}
    \begin{aligned}
        {\bf v} &= {\bf J}{\bf e} \\
        {\bf v} &= {\bf D} {\bf B}_\text{hb} \pamat {\bf e} \\
        \begin{bmatrix} v_q \\ v_p \end{bmatrix}
            &= \begin{bmatrix} D_{qq} & D_{qp} \\ D_{pq} & D_{pp} \end{bmatrix}
               \begin{bmatrix} B_{q\theta} & B_{q\phi} \\ B_{p\theta} & B_{p\phi} \end{bmatrix}
               \begin{bmatrix} P_{\theta x} & B_{\theta y} \\ B_{\phi x} & B_{\phi y} \end{bmatrix}
               \begin{bmatrix} e_x \\ e_y \end{bmatrix}
    \end{aligned}
\end{equation}

The \f${\bf D}\f$ matrix, if obtained using Hyperdrive, is already in the correct basis, and can be used as is.
However, if the RTS is used, the matrix must be permuted first.
In memory, matrices are represented as an array of \texttt{cuDoubleComplex}s (from the CUDA library), in the following order:
\begin{equation}
    \begin{bmatrix} 1 & 2 \\ 3 & 4 \end{bmatrix}.
\end{equation}
Therefore, the required permutation to get the {\bf D} matrix in the correct order is a ``reversal'' of the matrix elements in memory:
\begin{equation}
    \begin{aligned}
        \begin{bmatrix} 1 & 2 \\ 3 & 4 \end{bmatrix}
            &\rightarrow
            \begin{bmatrix} 4 & 3 \\ 2 & 1 \end{bmatrix} \\
        \begin{bmatrix} D_{pp} & D_{pq} \\ D_{qp} & D_{qq} \end{bmatrix}
            &\rightarrow
            \begin{bmatrix} D_{qq} & D_{qp} \\ D_{pq} & D_{pp} \end{bmatrix}
    \end{aligned}
\end{equation}

Currently, only Hyperbeam is used for the FEE beam model, and so the \f${\bf B}_{hp}\f$ is already given in the correct basis, \f$(\theta,\phi)\rightarrow(q,p)\f$.

Finally, the parallactic angle correction matrix is calculated by \vcsbeam{} in the \f$(x,y)\rightarrow(\theta,\phi)\f$ basis, as desired.

\subsection{Phasing up}

If the channelisation is fine enough, and if the tile separations (baselines) are not too long, then a time shift can be accurately effected by applying a phase ramp across frequency.

\subsection{Summing the voltages}

If we use subscripts \f$a\f$ and \f$f\f$ to represent tiles (``antennas'') and frequencies respectively, then the complete beamforming operation is
\begin{equation}
    {\bf e}_f = \sum_a e^{i\varphi_{a,f}} {\bf J}_{a,f}^{-1} {\bf v}_{a,f},
\end{equation}
where \f$\varphi_f\f$ is the frequency-dependent phase difference between the arrival of a signal at a given tile relative to an arbitrary reference point, or ``array phase centre''.

\section{Applications}

\section{Utilities}

\chapter{Implementation}

\section{Code glossary}

\subsection{\texttt{az}}
\begin{description}
    \item[Name:] \hyperlink{sec:coordslocalsky}{Geographic azimuth}\index{azimuth!geographic}
    \item[Units:] radians
    \item[Functions:] \hyperlink{fcn:calc_beam_geom}{\texttt{calc\_beam\_geom}}
\end{description}

\subsection{\texttt{el}}
\begin{description}
    \item[Name:] \hyperlink{sec:coordslocalsky}{Elevation}\index{elevation}
    \item[Units:] radians
    \item[Functions:] \hyperlink{fcn:calc_beam_geom}{\texttt{calc\_beam\_geom}}
\end{description}

\subsection{\texttt{za}}
\begin{description}
    \item[Name:] \hyperlink{sec:coordslocalsky}{Zenith angle}\index{zenith angle}
    \item[Units:] radians
    \item[Functions:] \hyperlink{fcn:calc_beam_geom}{\texttt{calc\_beam\_geom}}
\end{description}

\section{Function reference}

\subsection{\texttt{calc\_beam\_geom}}
\label{fcn:calc_beam_geom}

\appendix

\chapter{Description of RTS DI\_JonesMatrix file format}

\begin{lstlisting}
Franz Kirsten <franz.kirsten@curtin.edu.au>	17 June 2016 at 17:42
To: Sammy McSweeney <sammy.mcsweeney@gmail.com>, Steven Tremblay <Steven.Tremblay@curtin.edu.au>, Dr Stephen ORD <steve.ord@gmail.com>
Hi Sam,

Mitch finally got back to me on my questions about the individual entries in the DI_Jones matrices as output by the RTS. If you have a chance, could you go through those and try to get Andre's solutions to work once more? Are Mitch's answers enough or do you need more details?

It would be great if you could tell people next week that the beamformer also works with those tools now -- no pressure ;)

Cheers,
Franz




Hi Franz,

Sorry for the delay. Once you have this working I will let you know how to include the bandpass data.

Mitch

In an effort to automate calibrating the MWA for the tied-array beam, we
are currently trying to compare the the DI-Jones matrices from the RTS
with those produced by Andre Offringa's 'calibrate'. We know exactly
what is in Andre's solutions but are unsure about the
structur/formatting of the RTS solutions. Could you give a detailed
breakdown, please? In particular we were wondering about the following:

1) What is the number in the very first line?

This is the flux density of the calibrator that was used. You can ignore it.

The second line contains the model primary beam Jones matrix (in the direction of the calibrator).

2) Which column contains real/imaginary xx/yy/xy/yx?

XX_re XX_im XY_re XY_im YX_re YX_im YY_re YY_im

3) Are the numbers normalised in some way?

By the model primary beam Jones matrix (from the second line). If that is called B, and the direction-independent gain for tile i is G_i, then the numbers printed are G_i.B (call this J_i). So to get the direction-independent gain to compare with Andre's: G_i = J_i.inv(B).

4) Is there some conjugation going on?

No.
\end{lstlisting}

\chapter{Description of RTS BandpassCalibration file format}

\begin{lstlisting}
Sammy McSweeney <sammy.mcsweeney@gmail.com>	21 March 2017 at 12:09
To: Daniel.Mitchell@csiro.au, Franz Kirsten <franz.kirsten@curtin.edu.au>, Steven Tremblay <Steven.Tremblay@curtin.edu.au>
Hi Mitch,

Thanks for your explanation (via Franz, last year) of how the DIJones RTS output files are organised. If I could trouble you now for a similar explanation for how the Bandpass files are organised...
As far as I can make out, the format is something like:

freq00 freq01 freq02 ... (minus any flagged channels)
tilenum amp,phase, amp,phase, amp,phase, ...
tilenum amp,phase, amp,phase, amp,phase, ...
etc ...

I see that there are 8 rows per tile, but I'm not how those 8 break down (into polarisation combinations, for example). Incidentally, what do the "P" and "Q" mean in "PX", etc? (I see it cropping up in some of Steve Ord's plotting scripts.)

And, calling in your promised favour, anything you can tell me about "how to include the bandpass data" would be much appreciated.

Cheers,
~/Sam
\end{lstlisting}

\begin{lstlisting}
Daniel.Mitchell@csiro.au <Daniel.Mitchell@csiro.au>	23 March 2017 at 06:55
To: sammy.mcsweeney@gmail.com
Cc: franz.kirsten@curtin.edu.au, Steven.Tremblay@curtin.edu.au
Hi Sammy,

There are 8 lines for each antenna: Jm[0],Jf[0],Jm[1],Jf[1],Jm[2],Jf[2],Jm[3],Jf[3], where Jm contains the measured Jones matrices (measured separately for each frequency channel), and Jf contains low-order fits to these data, which are used in the RTS as the bandpass solutions (polynomials are fitted separately to the real and imaginary components of each Jones matrix element).

P and Q are the labels that we use for the two instrument polarisations (North-South and East-West respectively). X and Y are the standard sky polarisation coordinates, with X aligned with declination and Y aligned with RA.

The overall Jones matrix for a given channel of a given antenna will be the full-band matrix from the DIJones file multiplied by the bandpass matrix on the righthand side.

Cheers,
Mitch
\end{lstlisting}

\chapter{Description of Offringa/Hyperdrive calibration file format}
\label{app:offringa}

\begin{lstlisting}
Offringa's own documentation for the format of these binary files (copied from an email
dated 4 May 2016 from franz.kirsten@curtin.edu.au to sammy.mcsweeney@gmail.com). This is
itself copied from Andre Offringa's original C++ code, in his Anoko repository, in
mwa-reduce/solutionfile.h:

The solution file is used for storing calibration Jones matrices.
The format is as follows:
  Bytes |  Description
 -------+---------------
  0- 7  |  string intro ; 8-byte null terminated string "MWAOCAL"
  8-11  |  int fileType ; always 0, reserved for indicating something other than complex Jones solutions
 12-15  |  int structureType ; always 0, reserved for indicating different ordering
 16-19  |  int intervalCount ; Number of solution intervals in file
 20-23  |  int antennaCount ; Number of antennas that were in the measurement set (but were not necessary all solved for)
 24-27  |  int channelCount ; Number of channels in the measurement set
 28-31  |  int polarizationCount ; Number of polarizations solved for -- always four.
 32-39  |  double startTime ; Start time of solutions (AIPS time)
 40-47  |  double endTime ; End time of solutions (AIPS time)
 -------+-------------------
After the header follow 2 x nSolution doubles, with

nSolutions = nIntervals * nAntennas * nChannels * nPols

Ordered in the way as given, so:
double 0 : real of interval 0, antenna 0, channel 0, pol 0
double 1 : imaginary of interval 0, antenna 0, channel 0, pol 0
double 2 : real of interval 0, antenna 0, channel 0, pol 1
...
double 8 : real of interval 0, antenna 0, channel 1, pol 0
double nChannel x 8 : real of interval 0, antenna 1, channel 0, pol 0
etc.

ints are here always 32 bits unsigned integers, doubles are IEEE
double precision 64 bit floating points.
If a solution is not available, either because its data no data were
selected during calibration for this interval
or because the calibration diverged, a "NaN" will be stored in the
doubles belonging to that solution.
\end{lstlisting}

\chapter{Comparison of notation in other documents}

\section{Coordinate systems}

\begin{table}[!hb]
    \centering
    \caption{Comparison of coordinate notations}
    \label{tbl:notations}
    \begin{tabular}{l|cc|cc|cc}
        This document & \f$p\f$ & \f$q\f$ & \f$\theta\f$ & \f$\phi\f$ & \f$x\f$ & \f$y\f$ \\
        \hline
        MWA metafits files       & Y & X & - & - & - & - \\
        \citet{Sokolowski2017} & \f$y\f$ & \f$x\f$ & \f$\theta\f$ & \f$\phi\f$ & - & - \\
    \end{tabular}
\end{table}

\section{Jones matrices}
\label{sec:jones_matrices}

The notation used in this documentation most closely follows \citet{Ord2019}.

\begin{table}[!hb]
    \centering
    \caption{Comparison of matrix names}
    \label{tbl:notations}
    \begin{tabular}{l|ccc}
        This document & \f${\bf J}\f$ & \f${\bf D}\f$ & \f${\bf B}\f$ \\
        \hline
        \citet{Sokolowski2017} & \f${\bf J}\f$ & \f${\bf G}\f$ & \f${\bf E}\f$ \\
        \citet{Ord2019} & \f${\bf J}\f$ & \f${\bf D}\f$ & \f${\bf B}_s\f$/\f${\bf A}\f$ \\
    \end{tabular}
\end{table}

\printindex

\bibliography{vcsbeam_documentation}

\end{document}
