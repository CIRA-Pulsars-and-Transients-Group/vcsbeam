# Applications -- Make MWA Tied-Array Beam {#applicationsmakemwatiedarraybeam}

[TOC]

## Usage

usage: `make_mwa_tied_array_beam [OPTIONS]`

### Required options

| Short option | Long option | Description |
| --- | --- | --- |
| -c | --cal-metafits=FILE | Metafits file of the calibration observation |
| -C | --cal-location=PATH | The calibration solution file(s) |
| -m | --metafits=FILE | Metafits file of the target observation |
| -P | --pointings=FILE | File of pointings (RA and Dec) |

- For `-C`, the default calibration solution is a native mwa_hyperdrive FITS
  file. Hyperdrive writes one file per contiguous subband, so several files may
  be given (e.g. `-C solutions_*.fits`); each MPI task uses the file that covers
  its coarse channel. With `-O`, `-C` should instead be the full path to the
  single Offringa binary solution file.
- For `-P`, each line of the pointings file has the format
  `HH:MM:SS.S DD:MM:SS.S`, and the file may contain multiple pointings.

### Channelisation options

| Short option | Long option | Description | Default value |
| --- | --- | --- | --- |
| -A | --analysis_filter=FILTER | Fine channelisation filter | FINEPFB |
| -s | --smart | Use legacy fine channelisation | [off] |
| -S | --synth_filter=FILTER | Synthesis filter | LSQ12 |

- `-A` applies the named filter during fine channelisation (for MWAX only).
- `-s` uses the legacy settings for fine channelisation.
- `-S` applies the named filter during high-time resolution synthesis. FILTER
  can be MIRROR or LSQ12.

### Input options

| Short option | Long option | Description | Default value |
| --- | --- | --- | --- |
| -b | --begin=GPSTIME | Begin time (GPS seconds) | +0 |
| -d | --data-location=PATH | Input data directory | [current directory] |
| -f | --coarse-chan=CHAN | Coarse channel number | +0 |
| -T | --nseconds=VAL | Seconds of data to process | [as many as possible] |

- `-b` is the begin time of the observation, in GPS seconds. If GPSTIME starts
  with a '+' or a '-', then the time is taken relative to the start or end of
  the observation respectively.
- `-d` is the directory containing the recombined data.
- `-f` is the coarse channel number. If CHAN starts with a '+' or a '-', then
  the channel is taken relative to the first or last channel in the observation
  respectively. Otherwise, it is treated as a receiver channel number (0-255).
- `-T` processes VAL seconds of data.

### Output options

| Short option | Long option | Description | Default value |
| --- | --- | --- | --- |
| -p | --out-fine | Fine-channelised, full-Stokes PSRFITS output | [off] |
| -t | --max_t | Maximum seconds per output FITS file | 200 |
| -v | --out-coarse | Coarse-channelised, 2-pol (XY) VDIF output | [off] |

If neither `-p` nor `-v` is used, the default behaviour is to match the
channelisation of the input.

### Calibration options

| Short option | Long option | Description | Default value |
| --- | --- | --- | --- |
| -F | --flagged-tiles=FILE | FILE lists extra TileNames to flag | [none] |
| -O | --offringa | Solution is in Offringa binary format | [off] |
| -R | --ref-ant=TILENAME | Reference tile for phase rotation | [from file] |
| -U | --PQ-phase=PH,OFFS | Override the PQ phase correction | [from file] |
| -X | --cross-terms | Retain the PQ and QP terms | [off] |

- By default, tiles flagged in both the calibration and the observation
  metafits file are flagged in the beamformer; `-F` adds to these.
- `-O` selects the Offringa binary format instead of the default native
  mwa_hyperdrive FITS format.
- `-R` overrides the reference tile given in `pq_phase_correction.txt` for
  rotating the phases of the PP and QQ elements of the calibration solution. To
  turn off phase rotation altogether, set TILENAME=NONE.
- `-U` overrides the phase correction given in `pq_phase_correction.txt`. PH is
  given in rad/Hz and OFFS given in rad, such that the QQ element of the
  calibration Jones matrix for frequency F (in Hz) is multiplied by
  exp(PH\*F + OFFS). Setting PH = OFFS = 0 is equivalent to not performing any
  phase correction.
- `-X` retains the PQ and QP terms of the calibration solution.
- "[from file]" means the value in `pq_phase_correction.txt` is used.

### Memory options

| Short option | Long option | Description | Default value |
| --- | --- | --- | --- |
| -n | --nchunks=VAL | Processing chunks per second of data | 1 |

`-n` splits each second's worth of data into VAL processing chunks.

### Other options

| Short option | Long option | Description |
| --- | --- | --- |
| -h | --help | Print this help and exit |
| -V | --version | Print version number and exit |
