# Applications -- Make MWA Tied-Array Beam {#applicationsmakemwatiedarraybeam}

[TOC]

## Usage

usage: `make_mwa_tied_array_beam [OPTIONS]`

### Required options

| Short option | Long option | Description |
| ------------ | ----------- | ----------- |
| -c | --cal-metafits=FILE |  FILE is the metafits file pertaining to the calibration solution |
| -C | --cal-location=PATH |  PATH is the directory (RTS) or the file (OFFRINGA) containing the calibration solution |
| -m | --metafits=FILE     |  FILE is the metafits file for the target observation |
| -P | --pointings=FILE    |  FILE containing RA and Decs of multiple pointings in the format "HH:MM:SS.S DD:MM:SS.S ..." |

### Channelisation options

| Short option | Long option | Description | Default value |
| ------------ | ----------- | ----------- | ------------- |
| -A | --analysis_filter=FILTER | Apply the named filter during fine channelisation (for MWAX only) | FINEPFB |
| -s | --smart                  | Use legacy settings for fine channelisation | [off] |
| -S | --synth_filter=FILTER    | Apply the named filter during high-time resolution synthesis. FILTER can be MIRROR or LSQ12. | LSQ12 |

### Input options

| Short option | Long option | Description | Default value |
| ------------ | ----------- | ----------- | ------------- |
| -b | --begin=GPSTIME      | Begin time of observation, in GPS seconds.  If GPSTIME starts with a '+' or a '-', then the time is taken relative to the start or end of the observation respectively. | +0 |
| -d | --data-location=PATH | PATH is the directory containing the recombined data | [current directory] |
| -f | --coarse-chan=CHAN   | Coarse channel number. If CHAN starts with a '+' or a '-', then the channel is taken relative to the first or last channel in the observation respectively. Otherwise, it is treated as a receiver channel number (0-255) | +0 |
| -T | --nseconds=VAL       | Process VAL seconds of data | [as many as possible] |

### Output options

| Short option | Long option | Description | Default value |
| ------------ | ----------- | ----------- | ------------- |
| -p | --out-fine           |  Output fine-channelised, full-Stokes data (PSRFITS). If neither -p nor -v are used, default behaviour is to match channelisation of input. | [off] |
| -t | --max_t              |  Maximum number of seconds per output FITS file | 200 |
| -v | --out-coarse         |  Output coarse-channelised, 2-pol (XY) data (VDIF). If neither -p nor -v are used, default behaviour is to match channelisation of input. | [off] |

### Calibration options

| Short option | Long option | Description | Default value |
| ------------ | ----------- | ----------- | ------------- |
| -B | --bandpass           | Use the Bandpass (fine channel) as well as the DIJones (coarse channel) solutions (only relevant for RTS) | [off] |
| -F | --flagged-tiles=FILE | FILE is a text file including the TileNames of extra tiles to be flagged. By default, tiles flagged in both the calibration and the observation metafits file are flagged in the beamformer. The utility `rts_flag_ant_to_tilenames.py` can be used to convert the antenna numbers listed in the RTS `flagged_tiles.txt` file into human-readable tile names. | [no extra tiles are flagged] |
| -O | --offringa           | The calibration solution is in the Offringa format instead of the default RTS format. In this case, the argument to `-C` should be the full path to the binary solution file. | [off] |
| -R | --ref-ant=TILENAME   | Override the reference tile given in `pq_phase_correction.txt` for rotating the phases of the PP and QQ elements of the calibration solution. To turn off phase rotation altogether, set TILENAME=NONE. | [Use the value in `pq_phase_correction.txt`] |
| -U | --PQ-phase=PH,OFFS   | Override the phase correction given in `pq_phase_correction.txt`. PH is given in rad/Hz and OFFS given in rad, such that, the QQ element of the calibration Jones matrix for frequency F (in Hz) is multiplied by exp(PH\*F + OFFS). Setting PH = OFFS = 0 is equivalent to not performing any phase correction | [Use the value in `pq_phase_correction.txt`] |
| -X | --cross-terms        | Retain the PQ and QP terms of the calibration solution | [off] |

### Memory options

| Short option | Long option | Description | Default value |
| ------------ | ----------- | ----------- | ------------- |
| -n | --nchunks=VAL | Split each second's worth of data into VAL processing chunks | 1 |

### Other options

| Short option | Long option | Description |
| ------------ | ----------- | ----------- |
| -h | --help    | Print this help and exit |
| -V | --version | Print version number and exit |
