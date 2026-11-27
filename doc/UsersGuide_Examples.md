# User's Guide -- Examples {#usersguideexamples}

[TOC]

The example scripts use `module load vcsbeam/default`, which is usually set to
a sensible version on systems where VCSBeam is installed. If it is not, or if
you need a specific release, change it to whichever version is installed (see
`module avail vcsbeam`), for example `module load vcsbeam/v5.7`.

## Tied-array beamforming

### All 24 coarse channels

The script below forms a tied-array beam on all 24 coarse channels of an
observation, using one MPI task (and one GPU) per coarse channel.
Fill in the required inputs at the top, then submit it with `sbatch`.

**vcsbeam.sbatch**:

```
#!/bin/bash -l

#SBATCH --account=mwavcs-gpu
#SBATCH --partition=mwa-gpu
#SBATCH --job-name=beamform
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=24
#SBATCH --gpus-per-task=1
#SBATCH --time=04:00:00

module load vcsbeam/default

make_mwa_tied_array_beam -V

#===============================================================================
# Required inputs
#-------------------------------------------------------------------------------
# path to VCS metafits file
metafits=
# path to combined (.dat) or MWAX (.sub) data directory
datadir=
# path to calibration solution from hyperdrive
calsol=
# path to the calibrator observation metafits
calmetafits=
# the starting GPS second of the observation
startgps=
# how many seconds to process
duration=
# lowest coarse channel index
lowchan=
# path to pointings file
pointings=
#===============================================================================

srun make_mwa_tied_array_beam \
    -m ${metafits} \
    -b ${startgps} \
    -T ${duration} \
    -f ${lowchan} \
    -d ${datadir} \
    -P ${pointings} \
    -c ${calmetafits} \
    -C ${calsol} \
    -p -R NONE -U 0,0 -X -s
```

Notes on the script:

- `--ntasks=24` and `--gpus-per-task=1` give every coarse channel its own MPI
  task and GPU. Task `n` processes coarse channel `lowchan + n`, so the number
  of tasks sets the number of coarse channels that are beamformed.
- `calsol` is the native Hyperdrive FITS solution. Hyperdrive writes one file
  per contiguous subband, so `calsol` may also name several files (for example,
  `calsol="solutions_*.fits"`, left unquoted in the `srun` line so the shell
  expands it). Each task uses the file that covers its coarse channel.
- `-p` writes PSRFITS output (10 kHz channels, full Stokes).
- `-R NONE` and `-U 0,0` turn off the reference-tile phase rotation and the
  extra phase ramp, and `-X` keeps the cross terms of the calibration solution.
- `-s` turns on the SMART settings for the fine PFB.
- `datadir` can hold either combined (`.dat`) Legacy files or MWAX `.sub` files.
  The channel number in the `.sub` file names may or may not be zero-padded.
- Extra tiles can be flagged with `-F`, which takes a text file of tile names
  (one per line).
- The pointings file lists one pointing per line, as `hh:mm:ss.s dd:mm:ss.s`.

### A single second for testing

To test a setup quickly, use the same script with `duration=1`, and set
`startgps` to a second that contains a bright single pulse.

## Other processes

### Offline Fine PFB

The offline fine PFB creates recombined Legacy VCS files from MWAX data.
It processes one coarse channel at a time, so it is submitted as a job array.

**make_legacy.sbatch**:

```
#!/bin/bash -l

#SBATCH --account=mwavcs-gpu
#SBATCH --partition=mwa-gpu
#SBATCH --job-name=finepfb
#SBATCH --output=%x-%A_%a.out
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --array=109-132
#SBATCH --time=01:00:00

module load vcsbeam/default

# path to the MWAX data directory (and to the observation metafits)
datadir=
# the observation ID, which is also the starting GPS second
obsid=

srun fine_pfb_offline \
    -m ${datadir}/${obsid}.metafits \
    -b ${obsid} \
    -T 600 \
    -f ${SLURM_ARRAY_TASK_ID} \
    -d ${datadir}
```

### Legacy Offline Correlator

(This section is stale and the offline correlator is not actively used. The
module name, job settings and paths below may no longer be correct, and the
example is kept for reference only.)

The legacy offline correlator can be run on recombined files created with
Offline Fine PFB.

```
$ cd /astro/mwavcs/vcs/1320499816/cal/1320499816/vis
$ for freq in {109..132}
do
    sed "s/FREQ/${freq}/" < offline_correlate.sbatch > offline_correlate_${freq}.sbatch
    sbatch offline_correlate_${freq}.sbatch
done
```

**offline_correlate.sbatch**:

```
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mem=370gb
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --account=mwavcs
#SBATCH --export=NONE

module use /pawsey/mwa/software/python3/modulefiles
module load offline_correlator

for t in {8..599}
do

    gpssecond=$(echo "$t + 1320499816" | bc)
    gpubox=$(echo "FREQ - 108" | bc)

    INPUT_DATA_FILE=/astro/mwavcs/vcs/1320499816/combined/1320499816_${gpssecond}_chFREQ.dat
    START_SECOND=${gpssecond}
    DUMPS_PER_SECOND=1 # This sets the output time resolution
                    # (e.g. 20 --> 1/20 = 0.05s = 50 ms)
                    # Minimum allowed resolution is 2 ms
    CHANS_TO_AVERAGE=1 # This sets the output frequency resolution
                   # (e.g. 1 --> 1x10 kHz = 10 kHz)
    # This should be chosen to "match" the input channel.
    # This is not easy! (mwalib handles this, but at the
    # moment, offline_correlator is not using mwalib)
    GPUBOX_CHANNEL_NUMBER=${gpubox}
    OUTPUT_PREFIX=1320499816 # Output files begin with this

    srun -N 1 -n 1 offline_correlator \
        -d ${INPUT_DATA_FILE} \
        -s ${START_SECOND} \
        -r ${DUMPS_PER_SECOND} \
        -n ${CHANS_TO_AVERAGE} \
        -c ${GPUBOX_CHANNEL_NUMBER} \
        -o ${OUTPUT_PREFIX}

done
```

### Incoherent Beam

The incoherent beamformer has been tested and confirmed to work for both Legacy
and MWAX VCS data since the v5.6 release.

**vcsbeam_incoh.sbatch**:

```
#!/bin/bash -l

#SBATCH --account=mwavcs-gpu
#SBATCH --partition=mwa-gpu
#SBATCH --job-name=incoh
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=24
#SBATCH --gpus-per-task=1
#SBATCH --time=00:30:00

module load vcsbeam/default

# path to VCS metafits file
metafits=
# path to combined (.dat) or MWAX (.sub) data directory
datadir=
# the starting GPS second and how many seconds to process
startgps=
duration=
# lowest coarse channel index
lowchan=

srun make_mwa_incoh_beam \
    -m ${metafits} \
    -b ${startgps} \
    -T ${duration} \
    -f ${lowchan} \
    -d ${datadir} \
    -o incoh
```
