# CHANGELOG

<!-- markdownlint-disable MD024 -->
<!-- The same headings repeat for each release. -->

## v5.7

### Added

- MWAX `.sub` voltage files are found whether or not their channel number is
  zero-padded to three digits, so soft links are no longer needed for channels <
  100
- Native mwa_hyperdrive FITS calibration solutions can be read directly (now the
  default), without converting to the Offringa format. One file per subband may
  be given (e.g. `-C solutions_*.fits`)

### Changed

- Native mwa_hyperdrive FITS is now the default calibration solution format. The
  Offringa binary format is still available with `-O`
- Calibration messages from the Offringa and hyperdrive readers are now written
  with the MPI-aware logger, so they carry timestamps and MPI ranks. The
  solution file name is logged on its own line
- The per-second log message now reports the position of the coarse channel out
  of the total number of coarse channels (the MPI size), e.g. `[16/24]`, rather
  than `[N/1]`
- When MWAX `.sub` files are found only with unpadded channel numbers, the
  number of such files and (up to five of) their names are logged

### Removed

- Support for RTS calibration solutions (including the `-B/--bandpass` option,
  the `CAL_RTS` type and the `rts_flag_ant_to_tilenames.py` utility)
- The deprecated `mwa_plot_calibration` app

## v5.6 (2026-03-18)

### Added

- `-Z/--seconds-buffer=N` option to preload N seconds of data into memory (with
  OpenMP-parallel reads of multiple seconds) to improve I/O performance

### Changed

- Refactoring of the seconds-buffer implementation (now a named struct) and
  clearer comments and help text

### Fixed

- Memory access issues (segfaults) in the incoherent beamformer, including
  checks before freeing the calibration memory

## v5.5 (2025-10-21)

### Added

- `-D/--ds-factor` option to downsample the PSRFITS output (the incoherent
  beamformer never downsamples)
- `-k/--picket-fence` option to interpret the calibration solution appropriately
  for picket-fence data

## v5.4.1 (2025-09-26)

### Fixed

- Calibration channel indexing: the calibration channel and the processing
  channel indices were not aligned, so the solutions for only one coarse channel
  were applied across the whole bandwidth. This reduced the S/N, and is now
  corrected. Picket-fence processing was likely broken by this fix; native
  reading of Hyperdrive solutions (where the channel metadata is available) was
  identified as the long-term solution

## v5.4 (2025-08-18)

### Changed

- Updated the build script and module versions, and enabled unsafe atomics
  (`-munsafe-atomics`) so that the hardware implementation of atomic operations
  is used

## v5.3 (2025-05-21)

### Fixed

- Reverted the NaN guards added in v5.2.1 (calibration Jones matrices, tile sums
  and primary beam), as part of a temporary fix for NaN-filled output

## v5.2.1 (2025-04-29)

### Fixed

- Added guards against NaNs in the calibration Jones matrices, the voltage tile
  sums and the primary beam calculation, and a minor bug fix to the parallactic
  angle correction (the NaN guards were later removed in v5.3)

## v5.2 (2025-02-11)

### Fixed

- VDIF output: the real and imaginary parts were swapped

## v5.1 (2024-12-19)

### Changed

- Improved performance of the beamforming kernel: each beam is now computed by a
  single warp, with the sum performed as a parallel reduction across the warp,
  and just enough thread blocks are launched to avoid unnecessary context
  switches
- Removed a redundant input to `cu_invert_pfb`

### Fixed

- GPU macro definition for the device properties type

## v5.0.1-alpha (2024-11-19)

### Changed

- Updated the ROCm and hyperbeam API calls (hyperbeam > 0.7 and hipFFT for ROCm
  >= 5.7.3 are now required), and fixed the PAL CMake find recipe

### Fixed

- The Jones matrix configurations for each pointing in `vmCalcB` all pointed at
  the same memory, so all but one were in error

## v5.0.0-alpha (2024-09-03)

### Added

- GPU-agnostic code: VCSBeam can be built for both NVIDIA (CUDA) and AMD (HIP)
  GPUs (e.g. on Setonix), with `build_nvidia.sh` and `build_hip.sh` build
  scripts, a single refactored top-level `CMakeLists.txt`, and GPU macros in
  `gpu_macros.h`
- `-N/--out-nstokes` option to choose the number of Stokes parameters in the
  PSRFITS output (1 for Stokes I only, or 4 for IQUV)
- Command-line arguments and error checking in the RTS tile conversion script
- `nchunks` option to divide the data into chunks to fit onto smaller GPUs, and
  `NPOL` in the VDIF header

### Fixed

- Frequency table mislabelling in PSRFITS output for picket-fence data
- The `orig_nbits` PSRFITS header entry was not set, which caused psrfits_utils
  to write all zeroes
- Start time (including offsets) in the output file headers, for both VDIF and
  PSRFITS
- Reading of Offringa-style calibration solutions, including the number of
  channels used, and GPU memory address alignment errors
- VDIF output bugs (including which half of the buffer was filled) and the VDIF
  ASCII header, so that DSPSR can read it

## v4.2 (2023-09-09)

### Added

- Re-implemented multipixel beamforming ("Multibeam functionality"; the first
  non-crashing version, not yet tested on multiple pulsars)
- mwalib versioning information in the docs, and updated examples in the
  Beamforming and Examples user guides

### Changed

- Reworked the beamformer, calibration and Jones-matrix code (`form_beam.cu`,
  `calibration.c`, `jones.c`, `metadata.c`) and the `vcsbeam.h` API to support
  multipixel beamforming
- Some warning messages now go to stdout

### Fixed

- A channel-indexing bug when reading Offringa calibration solution files
- Changed `DAT_OFFS` to its minimum value (a possible bug fix, untested)

## v4.1 (2023-06-16)

### Added

- Error checking for the read-buffer states `VM_END_OF_DATA`,
  `VM_READ_BUFFER_NOT_SET` and `VM_READ_BUFFER_LOCKED`, via `vmCheckError`

### Changed

- Updated the calibration example images and text in the docs, and added mwalib
  versioning information to the Installation docs
- Removed some debugging output and fixed a documentation typo (#20)

### Fixed

- Non-deterministic beamformer output, caused by a race condition in the CUDA
  beamforming kernel while antennas are summed (`form_beam.cu`)
- Updated `mwa_track_primary_beam_response` and `mwa_tied_array_beam_psf` to
  match the memory requirements of the latest hyperbeam version

## v4.0 (2022-10-23)

### Added

- Usage documentation for several more apps, and PSRFITS output changes in
  `beam_psrfits.c`

### Changed

- Updated to hyperbeam v0.5.0 (including changes to `primary_beam.c`)
- The start frequency used for counting coarse and fine channels is now absolute
- Updated the build and find modules and the README (Codacy)

### Removed

- The `plot_calibration.py` utility and its documentation page, because
  hyperdrive plots calibration solutions better

### Fixed

- A renormalisation bug, so the output is now correct
- An incorrect line in the error handling in `primary_beam.c`
- An installation bug caused by the removed utility

## v3.4 (2022-08-08)

### Added

- Options in `mwa_tied_array_beam_psf` so that the image can be made anywhere in
  the sky, not just centred on the tied-array beam
- A large expansion of the Doxygen documentation: docstrings for all functions,
  a User's Guide, beamforming and calibration workflow pages with graphviz
  diagrams, an Offline correlator page skeleton, and a VCSBeam logo

### Changed

- Updated to the new mwa_hyperbeam API, and updated the docs for the hyperbeam
  version number

### Fixed

- Several `mwa_plot_calibration` bugs: calibration corrections are now applied
  to the printed output, it assumes 128 fine channels, RTS file reading is
  fixed, and a missing usage option is added

## v3.3 (2022-01-03)

### Added

- Doxygen documentation, with docstrings across the source files (`jones.c`,
  `performance.c`, `pfb.cu` and others). The existing high-level docs were
  converted to Doxygen/markup and material from the MWA telescope wiki was added
- A GitHub workflow that publishes the docs to github.io

### Removed

- Unused functions in `beam_vdif.c` and leftover backup files

### Fixed

- PSRFITS header frequency labelling (Hz changed to MHz)
- A segfault when a flagged tile name does not exist, and a crash when no
  flagged-tilenames file is given

## v3.2 (2021-12-21)

### Added

- Beamforming directly on MWAX data, without writing intermediate files first
  (first light on 2021-12-21)
- `-s` option in the tied-array beamformer to turn on SMART settings for the
  fine PFB
- Chunking in the fine PFB, and a generalised forward PFB that no longer needs
  to match the FPGA settings (with a scaling factor so that its output matches
  the FPGA version)
- Improved choice between PSRFITS and VDIF output in the tied-array beamformer
- Performance statistics now report only the stopwatches that were used, and are
  printed more clearly

### Changed

- Refactored the code into a `vcsbeam_context` API, with data buffer management,
  chunk-size calculations and datatypes moved into it. The old
  `vcsbeam_metadata` was renamed, and unneeded GPU mallocs per second of input
  were removed
- Removed `ANTENNA_ORDER.md` and fixed documentation typos

### Fixed

- Re/im ordering in the fine PFB output, including the full-precision mode
- PSRFITS header and frequency labelling bugs, and wrong channel counts in the
  MWAX case
- The reference tile phase correction is now applied to qq instead of pp
- Log output for non-MPI runs

## v3.1 (2021-12-03)

### Changed

- Fixed incorrect information in the documentation and added more text

### Fixed

- Restored the reference tile phase removal so that each element of the Jones
  matrix is treated independently, as before

## v3.0 (2021-12-01)

### Changed

- Rewrote the matrix basis transformations in the docs as (c,d)<--(a,b), and
  added documentation on Hyperdrive, the Offringa calibration solution file
  format and beamforming

### Fixed

- The polarisation and conjugation handling. The RTS does not conjugate the sky,
  because swapping the real and imaginary parts of the voltages is equivalent,
  so the "RTS conjugates the sky" workaround was removed from the beamformer
- The parallactic angle correction, so that the final coherency matrix has sky
  XX in the top-left corner
- Added a matrix column swap (`swapcols`) for Hyperdrive solutions and reduced
  the row swaps/reversals needed to make RTS and Offringa-style solutions
  consistent
- Known issue at this release: some polarisation leakage (too much Stokes V on
  some test observations) was still being investigated; Stokes I looks good

## v2.20 (2021-11-29)

### Added

- `plot_calibration.py` now uses argparse and writes PNG files instead of
  plotting to screen, and the calibration plots are much improved
- `mwa_plot_calibration` can print user-defined strings (e.g. `nan`) for zero
  matrices

### Changed

- Reworked the order in which RTS calibration solutions are constructed, which
  changes the output
- The "divide by reference tile" step now divides all four 2x2 Jones elements by
  the phase of just one of them
- Settled the X/Y polarisation labelling and notation in the documentation (for
  legacy data the metafits labels are not what they appear to be), and renamed
  some matrices in the docs

### Fixed

- Calibration solutions failing for the MWAX test observation (obsid 1320499816)
- Documentation fixes: calibration section text, a missing `.bib` file, a
  corrected equation and a code comment

## v2.19 (2021-11-24)

### Added

- Tile flagging during beamforming, using the flagged-tile information
- `-B` option to make inclusion of the bandpass solutions optional; the bandpass
  solutions are reinstated in the calibration
- Reading of Offringa-style calibration solutions, implemented and used in
  `mwa_plot_calibration`
- A more flexible `mwa_track_primary_beam_response`: output format, optional
  primary beam Jones matrices, customisable start/end times and time stride, and
  an option to turn the parallactic angle correction on or off
- New Python utilities for interactively finding the best phase slope
  (`pol_correction.py`) and for printing metafits inputs; utils are now
  installed via a CMake file
- The PSF app prints all Stokes responses
- A LaTeX/PDF documentation of the notation and calibration was started, and the
  primary beam matrix operations were corrected to match it

### Fixed

- Non-master MPI tasks did not wait for their data to be spliced
- `mwa_plot_calibration` used the wrong metafits file for the calibration
  correction
- A bandpass channel-counting bug, a single-antenna case in
  `rts_flag_ant_to_tilenames.py`, and hardcoded filenames in `pol_correction.py`

## v2.18 (2021-11-03)

### Added

- New app `mwa_tied_array_beam_psf`, which computes a thumbnail image of the
  tied-array beam PSF for an observation
- New app `mwa_plot_calibration` with the accompanying `plot_calibration.py`,
  which produce calibration plots (including an `-N` option)
- A utility for extracting metadata from metafits files
- Customisable calibration correction tools on the command line, and custom
  antenna flagging in the tied-array beamformer
- The tied-array beamformer prints the mwalib version in its output
- A documentation page on antenna ordering

### Changed

- The PQ phase correction is reworked, with an extra column in
  `pq_phase_correction.txt`, and is now actually applied. The phase slope is
  applied to PP instead of QQ
- The bandpass solutions were disabled again in this release

### Fixed

- The antenna/polarisation ordering problem (legacy order is Y,X rather than
  X,Y) that stopped the Vela observation being detected, and the order of the
  complex weights
- A segfault when there were fewer calibration channels than expected
- The antenna order of the offline fine PFB output

## v2.17 (2021-09-26)

### Changed

- Dependency handling is more flexible and robust (not yet tested on systems
  where some dependencies are missing)

## v2.16 (2021-09-25)

### Added

- `offline_correlator` as a separate directory with its own CMake build
- xGPU is now an optional compile-time dependency
- README updates for hyperbeam and VDIFIO

### Fixed

- Restored the runtime files that were lost in the move to a library
- Removed C99 non-compliant code from `offline_correlator` and a redundant `-i`
  option

## v2.15 (2021-09-22)

### Added

- The code is now built as a library (CMake generates `vcsbeam.h` from
  `vcsbeam.h.in`), with executables linked against it
- PFB flags let users manage their own memory allocation

### Changed

- VDIFIO is now a full dependency pulled from its GitHub source. The VDIF header
  fix is included
- Removed a duplicate allocation function from `form_beam.c`

## v2.14 (2021-09-16)

### Added

- Forward (fine) PFB on the GPU (weighted overlap-add), with a new
  `fine_pfb_offline` program that writes VCS-recombined-style files. It
  installs, and the output sample position was corrected
- A new `vcsbeam_metadata` struct that abstracts the mwalib metadata handling,
  now used in the tied-array and incoherent beamformers

### Changed

- `ipfb` is renamed to `pfb`, and filter coefficients are loaded via a new
  struct
- Antenna flagging now uses mwalib's built-in list
- `xy_phase_correction` is renamed `pq_phase_correction`
- Runtime files install to a different default directory
- The calibration `D` array is now 1D

## v2.13 (2021-09-06)

### Changed

- Tidied the tied-array beamformer options and code, and moved the write step
  before the read step to maximise MPI concurrency
- Removed the legacy `splice_psrfits` sources

### Fixed

- The remaining MPI splicing artefacts (rescaling of channels kernel)

## v2.12 (2021-09-06)

### Added

- MPI support in both `make_mwa_incoherent_beam` and `make_mwa_tied_array_beam`:
  channels are distributed over MPI ranks and the PSRFITS output is spliced via
  MPI (`beam_psrfits`). The MPI rank is added to the log output
- `mwa_track_primary_beam_response` gains the distance from the tile pointing
  centre and an optional array-factor calculation

### Changed

- Incoherent beamforming was removed from the tied-array beamformer. The number
  of coarse channels processed must equal the number of MPI jobs, and `-F` again
  means the actual number of channels
- Compatible with mwalib v0.10.0, and the CMake output is reordered and builds
  more portably (fewer `strdup` calls)
- The incoherent sum discretisation was merged into the flattened-bandpass
  kernel

### Fixed

- A race between MPI processes altering the data buffer, and an indexing bug
  that affected the last sample of every second in the incoherent beamformer

## v2.11 (2021-08-25)

### Added

- New utility `mwa_track_primary_beam_response`
- Both beamformers now loop over coarse channels (mwalib limits this to one in
  practice)

### Changed

- Large clean-up: `beam_common` retired and its functions re-homed, plus renamed
  files and macros. The `begin` option is no longer required

## v2.10 (2021-08-23)

### Added

- Option for a custom output filename in the incoherent beamformer
- The bandpass solutions are put back in and work with the PFB inversion

### Changed

- The number of seconds to process replaces the end GPS time on both command
  lines, and the `begin` option now uses a new style
- The `t` option becomes `T`, and the default data directory is the current
  directory
- Functions are agnostic to the number of pointings

### Fixed

- Tiles with more than two dead dipoles are now included in beamforming

## v2.9 (2021-08-21)

### Added

- Stand-alone incoherent beam program `make_mwa_incoh_beam` with a new CUDA
  incoherent beam kernel; `performance.c` holds its runtime calculations and
  logging
- `xy_phase_correction.txt` runtime file for the XY phase correction, installed
  automatically. If it is missing, the beamformer continues without the
  correction

### Changed

- `make_beam` is renamed `make_mwa_tied_array_beam`, and pointings must be
  supplied via a file
- The PFB inversion kernel is significantly faster
- Primary beam (`primary_beam.c`) and geometric delay (`geometric_delay.c`) code
  is extracted from `get_jones`, with other calibration code moved to
  `calibrate.c` and mwalib wrappers to `metadata.c`

## v2.8 (2021-08-18)

### Added

- Bandpass calibration solutions can be read from the calibration metafits file

### Changed

- Calibration-flagged tiles are added to the overall flagged tiles
- `get_delays` is renamed `get_jones` (`get_delays_small.c` becomes
  `get_jones.c`), and its geometric calculations are separated out
- The alignment matrix is applied in `calibrate.c`
- Removed the `-J` and `-B` options and fully retired the old `metafits_info`
  struct in favour of mwalib. The array dimensions in the CUDA code now come
  from mwalib constants

## v2.7 (2021-08-13)

### Added

- vcsbeam extracted from the vcstools repository as its own repo, with a README
  and LICENSE. FFTW and the deprecated PFB inversion functions were removed
- Switched to mwalib for metadata, antenna ordering, delays and amplitudes, and
  calibration
- Options in `make_beam` for applying a UV phase correction, and for John's
  method of removing calibration artefacts (default)
- `sky_model_jones_comparison` and a first draft of `azza_track` tools
- Command-line option for the calibration metafits file

### Changed

- Complex numbers converted to CUDA complex types and `mycomplex.h` removed.
  SLALIB replaced with PAL
- The analytic beam model was removed. The FEE2016 beam application was
  simplified, and its Jones frame matched to the analytic beam conventions
- The reference antenna division now divides all four Jones elements
- The calibration channel width is now read from the metafits

### Fixed

- A UTC time error in `process_vcs.py`, which put the LST an hour out, and a
  Jones matrix memory-free issue
- Handling of dipole configurations with three or more dead dipoles (untested)

## v2.6 (2021-04-29, also tagged v2.5)

### Fixed

- Incoherent beamforming now works with SSDs

## v2.4 (2020-11-06)

### Changed

- Switched the astronomy library calls in `get_delays_small.c`, `beam_psrfits.c`
  and `beam_vdif.c` from SLALIB to PAL, with matching `CMakeLists.txt` changes.
  The PAL port compiles, but its commits note that it has not been tested for
  correctness

### Fixed

- Two bugs introduced by merge-conflict resolution

## v2.3 (2020-11-06)

### Added

- Support for the MWA FEE 2016 beam model through `mwa_hyperbeam`, selected with
  the `-H` option (which is checked, and also accepts "analytic")
- `--analytic-beam` option (renamed from `--beam-model`)
- The FEE beam model is read once at start-up, which makes it much faster
- Dipole delays and amplitudes are now read per tile, so each tile gets an
  accurate Jones matrix
- Repeated dipole configurations (for example dead dipoles) are no longer
  recalculated, via a new `hash_dipole_config` function

### Changed

- Added a cleanup and comments for `calcEjones`, the analytic beam model in
  `get_delays_small.c`
- Debug output of tile delays was removed from the `make_beam` output

### Fixed

- A GPU "illegal memory access" error and a bug that could flag more tiles than
  necessary
- A type mismatch in the `metainfo` struct in `beam_common.c`
- FEE beams are now recalculated for each pointing in multipixel mode
- The X and Y definitions were swapped, which corrected the signs of the
  polarisation angle slope and Stokes V

## v2.2 (2020-09-02)

### Added

- `make_beam` now checks how much GPU memory is available and splits the
  processing into chunks that fit
- An option to set how much GPU memory the beamformer may use, with help text
  for the `-g` option
- An option to set the maximum seconds per output FITS file (default 200 s)

### Changed

- Reworked how timing (benchmarking) is handled in the beamformer
- Raised the initial limit used in the GPU memory check

## v2.1 (2020-06-24)

### Added

- `-S [synth_filter]` option in `make_beam` to choose the inverse PFB synthesis
  filter
- Least-squares (12-tap) synthesis filter coefficients in `filter.h`
- `-A` option, so that a single RF stream can be inverted without beamforming

### Changed

- Rewrote the inverse PFB code in `ipfb.cu` with a new synthesis filter kernel
- Raised the VDIF output scale factor from 10 to 15 to make better use of the
  available bits
- Updated the `make_beam` usage text

### Fixed

- A redundant free that caused segmentation faults
- Compilation errors and warnings

## v2.0 (2020-04-21)

### Changed

- Swapped the cross terms in the beamformer calibration step. This is the
  correct form, and gives a higher signal-to-noise ratio for detections

### Fixed

- The first 415 samples of each second being full of RFI: the bandpass was being
  flattened more than once

## v1.8 (2020-02-25, also tagged v1.7)

### Added

- An option to run the beamformer in serial or with OpenMP
- Mean and standard deviation in the section timings

### Changed

- Calibration is now applied once, in a split beamforming kernel, rather than
  for every pointing. This makes the GPU calculation about 63% faster after the
  first pointing
- Bandpass flattening now happens once instead of in chunks
- Removed the unnecessary pointing dimension from the Jones matrix
- Changed build flags, including dropping `-lineinfo`
- GPU work was found to be slower inside OpenMP loops, so it was made serial
  again

### Fixed

- A VDIF bug that slowly lowered the gain over long observations
- Unused-variable warnings in the inverse PFB, and the scope of `start` in
  `make_beam.c` (a Codacy fix)

## v1.6 (2019-10-16)

### Changed

- Updated the CMake files and the beamformer to handle the new version number
  format

## v1.5 (2019-09-27, also tagged v1.5.0)

### Added

- First tagged release of the MWA beamformer, with `make_beam` and `get_delays`.
  They read RTS and Offringa-style calibration solutions and write PSRFITS and
  VDIF output, with `-V` printing the version
- A fringe builder, options to swap or conjugate polarisation and sky, and a
  `-G` option
- Options to set the begin and end times (`-b`, `-e`) and the coarse channel
  (`-N`). Parallelism moved out of `make_beam` and into a wrapper script
- A heap-based mode as an alternative to temporary files, and an AGC rate option
- `antenna_map.txt` in place of the hard-coded antenna ordering, and custom tile
  flag files (also `-F` in `process_vcs.py`)
- A PSRFITS full-polarisation builder, simultaneous PSRFITS and VDIF output, and
  an option for summed-polarisation (Stokes I) output
- A multi-pixel beamformer with up to 120 pointings, each written to its own
  sub-directory, in both PSRFITS and VDIF
- Asynchronous read, calculate and write (about 1.7x faster), plus a README and
  AFL-3.0 licence information

### Changed

- Cleaned up the `CMakeLists.txt` files, standardised version handling, and
  removed the deprecated beamformer code (v1.4.0). The slimmed-down
  `make_beam_small` was renamed `make_beam`
- Moved bandpass flattening and the beamformer kernel to the GPU and optimised
  them. Memory is now allocated once
- Reworked `get_delays` so that it no longer writes a Jones file and normalises
  the Jones matrices differently, which means results differ from earlier
  versions
- Reduced pointing output precision to 0.01 arcsec

### Fixed

- NaN and zero output on Athena, by initialising GPU arrays to zero
- A polarisation bug where PSRFITS output reported one polarisation instead of
  four
- Several VDIF, normalisation, memory and pointer-arithmetic bugs
