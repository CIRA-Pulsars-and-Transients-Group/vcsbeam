# CHANGELOG

## v4.2

### Added

- Multibeam functionality
- MWAX `.sub` voltage files are found whether or not their channel number is zero-padded to three digits, so soft links are no longer needed for channels < 100
- Native mwa_hyperdrive FITS calibration solutions can be read directly (now the default), without converting to the Offringa format. One file per subband may be given (e.g. `-C solutions_*.fits`)

### Removed

- Support for RTS calibration solutions (and the `-B/--bandpass` option), the `rts_flag_ant_to_tilenames.py` utility, and the deprecated `mwa_plot_calibration` app
