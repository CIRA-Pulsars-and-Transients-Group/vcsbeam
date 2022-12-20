# Change Log

## [unreleased]

### Added

- [ ] Added writing out support for filterbank (`.fil`) files
- [ ] Added writing out support for HDF5 files
- [ ] Added downsampling in time and/or frequency before writing out
- [ ] Write out (many) coherent dedispersion trials

### Changed

- [x] `make_mwa_tied_array_beam` options `-p` (PSRFITS output) and `-v` (VDIF OUTPUT) replaced by more generalisable `-o` option.
- [x] Made `make_mwa_tied_array_beam`'s dependence on `PSRFITS_UTILS` optional.

## v4.0

### Fixed

- Fixed channel normalisation bug (thanks [Jared Moseley](https://github.com/Jared-Moseley)!)
