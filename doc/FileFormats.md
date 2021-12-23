# File Formats {#fileformats}

[TOC]

## Data files

### Legacy recombined files

@todo Write description of the legacy recombined file format

### MWAX files

A full description of the MWAX hight time resolution format can be found on [this MWA Telescope Wiki page](https://wiki.mwatelescope.org/display/MP/MWA+High+Time+Resolution+Voltage+Capture+System).

## Calibration solutions

### RTS DI\_JonesMatrix file format {#rtsfileformat}

```
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
```

### RTS BandpassCalibration file format

```
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
```

```
Daniel.Mitchell@csiro.au <Daniel.Mitchell@csiro.au>	23 March 2017 at 06:55
To: sammy.mcsweeney@gmail.com
Cc: franz.kirsten@curtin.edu.au, Steven.Tremblay@curtin.edu.au
Hi Sammy,

There are 8 lines for each antenna: Jm[0],Jf[0],Jm[1],Jf[1],Jm[2],Jf[2],Jm[3],Jf[3], where Jm contains the measured Jones matrices (measured separately for each frequency channel), and Jf contains low-order fits to these data, which are used in the RTS as the bandpass solutions (polynomials are fitted separately to the real and imaginary components of each Jones matrix element).

P and Q are the labels that we use for the two instrument polarisations (North-South and East-West respectively). X and Y are the standard sky polarisation coordinates, with X aligned with declination and Y aligned with RA.

The overall Jones matrix for a given channel of a given antenna will be the full-band matrix from the DIJones file multiplied by the bandpass matrix on the righthand side.

Cheers,
Mitch
```

### Offringa/Hyperdrive calibration file format {#offringa}

```
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
```
