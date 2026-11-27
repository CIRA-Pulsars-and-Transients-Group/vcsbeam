# CHANGELOG

## v5.7

### Added

- MWAX `.sub` voltage files are found whether or not their channel number is zero-padded to three digits, so soft links are no longer needed for channels < 100
- Native mwa_hyperdrive FITS calibration solutions can be read directly (now the default), without converting to the Offringa format. One file per subband may be given (e.g. `-C solutions_*.fits`)

### Changed

- Native mwa_hyperdrive FITS is now the default calibration solution format. The Offringa binary format is still available with `-O`
- Calibration messages from the Offringa and hyperdrive readers are now written with the MPI-aware logger, so they carry timestamps and MPI ranks. The solution file name is logged on its own line
- The per-second log message now reports the position of the coarse channel out of the total number of coarse channels (the MPI size), e.g. `[16/24]`, rather than `[N/1]`
- When MWAX `.sub` files are found only with unpadded channel numbers, the number of such files and (up to five of) their names are logged

### Removed

- Support for RTS calibration solutions (including the `-B/--bandpass` option, the `CAL_RTS` type and the `rts_flag_ant_to_tilenames.py` utility)
- The deprecated `mwa_plot_calibration` app

## v4.2

### Added

- Multibeam functionality
