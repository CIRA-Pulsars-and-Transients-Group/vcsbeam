# Applications -- Fine PFB Offline {#applicationsfinepfboffline}

[TOC]

## Usage

usage: `fine_pfb_offline [OPTIONS]`

### Required options

| Short option | Long option | Description |
| --- | --- | --- |
| -m | --metafits=FILE | Metafits file of the target observation |

### Optional options

| Short option | Long option | Description | Default value |
| --- | --- | --- | --- |
| -A | --analysis_filter=FILTER | Fine channelisation filter | FINEPFB |
| -b | --begin=GPSTIME | Begin time (GPS seconds) | +0 |
| -d | --data-location=PATH | Input data directory | [current directory] |
| -f | --coarse-chan=CHAN | Coarse channel number | +0 |
| -n | --nchunks=VAL | Processing chunks per second of data | 1 |
| -T | --nseconds=VAL | Seconds of data to process | [as many as possible] |

- `-A` applies the named filter during fine channelisation. The file
  `[RUNTIME_DIR]/FILTER.dat` must exist.
- `-b` is the begin time of the observation, in GPS seconds. If GPSTIME starts
  with a '+' or a '-', then the time is taken relative to the start or end of
  the observation respectively.
- `-d` is the directory containing the recombined data.
- `-f` is the coarse channel number. If CHAN starts with a '+' or a '-', then
  the channel is taken relative to the first or last channel in the observation
  respectively. Otherwise, it is treated as a receiver channel number (0-255).
- `-n` splits each second's worth of data into VAL processing chunks.
- `-T` processes VAL seconds of data.

### Other options

| Short option | Long option | Description |
| --- | --- | --- |
| -h | --help | Print this help and exit |
| -V | --version | Print version number and exit |
