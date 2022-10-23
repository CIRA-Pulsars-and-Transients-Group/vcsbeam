# Applications -- Fine PFB Offline {#applicationsfinepfboffline}

[TOC]

## Usage

usage: `fine_pfb_offline [OPTIONS]`

### Required options

| Short option | Long option | Description |
| ------------ | ----------- | ----------- |
| -m | --metafits=FILE     |  FILE is the metafits file for the target observation |

### Optional options

| Short option | Long option | Description | Default value |
| ------------ | ----------- | ----------- | ------------- |
| -A | --analysis_filter=FILTER | Apply the named filter during fine channelisation. File [RUNTIME_DIR]/FILTER.dat must exist | FINEPFB |
| -b | --begin=GPSTIME      | Begin time of observation, in GPS seconds.  If GPSTIME starts with a '+' or a '-', then the time is taken relative to the start or end of the observation respectively. | +0 |
| -d | --data-location=PATH | PATH is the directory containing the recombined data | [current directory] |
| -f | --coarse-chan=CHAN   | Coarse channel number. If CHAN starts with a '+' or a '-', then the channel is taken relative to the first or last channel in the observation respectively. Otherwise, it is treated as a receiver channel number (0-255) | +0 |
| -n | --nchunks=VAL | Split each second's worth of data into VAL processing chunks | 1 |
| -T | --nseconds=VAL       | Process VAL seconds of data | [as many as possible] |

### Other options

| Short option | Long option | Description |
| ------------ | ----------- | ----------- |
| -h | --help    | Print this help and exit |
| -V | --version | Print version number and exit |
