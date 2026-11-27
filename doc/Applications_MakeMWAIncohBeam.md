# Applications -- Make MWA Incoh Beam {#applicationsmakemwaincohbeam}

[TOC]

## Usage

usage: `make_mwa_incoh_beam [OPTIONS]`

### Required options

| Short option | Long option | Description |
| --- | --- | --- |
| -m | --metafits=FILE | Metafits file of the target observation |

### Input & output options

| Short option | Long option | Description | Default value |
| --- | --- | --- | --- |
| -b | --begin=GPSTIME | Begin time (GPS seconds) | +0 |
| -d | --data-location=PATH | Input data directory | [current directory] |
| -f | --coarse-chan=CHAN | Coarse channel number | +0 |
| -o | --outfile | Base name of the output PSRFITS file | see below |
| -S | --max_output_t=SECS | Maximum seconds per output FITS file | 200 |
| -T | --nseconds=VAL | Seconds of data to process | [as many as possible] |

- `-b` is the begin time of the observation, in GPS seconds. If GPSTIME starts
  with a '+' or a '-', then the time is taken relative to the start or end of
  the observation respectively.
- `-d` is the directory containing the recombined data.
- `-f` is the coarse channel number. If CHAN starts with a '+' or a '-', then
  the channel is taken relative to the first or last channel in the observation
  respectively. Otherwise, it is treated as a receiver channel number (0-255).
- `-o` is the base name for the output PSRFITS file. The default is
  `<PROJECT>_<OBSID>_incoh_ch<CHAN>`.
- `-S` is the maximum number of SECS per output FITS file.
- `-T` processes VAL seconds of data.

### Other options

| Short option | Long option | Description |
| --- | --- | --- |
| -h | --help | Print this help and exit |
| -V | --version | Print version number and exit |
