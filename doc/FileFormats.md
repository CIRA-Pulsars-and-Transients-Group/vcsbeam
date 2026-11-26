# File Formats {#fileformats}

[TOC]

## Data files

### Legacy recombined files

@todo Write description of the legacy recombined file format

### MWAX files

A full description of the MWAX hight time resolution format can be found on [this MWA Telescope Wiki page](https://wiki.mwatelescope.org/display/MP/MWA+High+Time+Resolution+Voltage+Capture+System).

## Calibration solutions

### Hyperdrive FITS calibration file format {#hyperdrivefileformat}

Native `mwa_hyperdrive` solutions (the default) are FITS files read directly by VCSBeam. The relevant HDUs are:

| HDU | Contents |
| --- | --- |
| `SOLUTIONS` | Double-precision image with FITS axes `8 x nchan x ntiles x ntimeblocks`. The 8 values are the real and imaginary parts of the Jones matrix elements (PP, PQ, QP, QQ). NaN means the solution is unavailable. Only the first timeblock is used. |
| `TILES` | Binary table. The `Antenna` column (metafits antenna number) and `Flag` column (non-zero means flagged) are used. |
| `CHANBLOCKS` | Binary table. The `Freq` column (channel centre, Hz) is used to find the channels belonging to each coarse channel, and `Flag` marks flagged channels. |

Hyperdrive writes one file per contiguous subband. `-C` can be followed by several files (e.g. a shell glob, `-C solutions_*.fits`); each MPI task uses the first file whose frequencies cover its coarse channel, so any mix of subband widths (up to one file per coarse channel for picket-fence data) is supported.

### Offringa/Hyperdrive calibration file format {#offringafileformat}

Offringa's own documentation for the format of these binary files (copied from an email dated 4 May 2016 from `franz.kirsten@curtin.edu.au` to `sammy.mcsweeney@gmail.com)`.
This is itself copied from Andre Offringa's original C++ code, in his Anoko repository, in `mwa-reduce/solutionfile.h`:
```
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
```
