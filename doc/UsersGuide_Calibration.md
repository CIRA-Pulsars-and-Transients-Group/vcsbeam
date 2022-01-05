# User's Guide -- Preparing a calibration solution {#usersguidecalibration}

[TOC]

## RTS vs Hyperdrive

Which calibration software to use (RTS vs Hyperdrive) depends on whether the data is Legacy or MWAX, and whether it is VCS or already-correlated data (e.g. a dedicated calibration observation).
The following table summarises the possibilities:

[RTS]: https://wiki.mwatelescope.org/display/MP/Documentation#Documentation-CalibratingwiththeRealTimeSystem(RTS)
[Hyperdrive]: https://wiki.mwatelescope.org/pages/viewpage.action?pageId=52068764

|        | VCS                      | Correlated visibilities              |
| ------ | ------------------------ | ------------------------------------ |
| Legacy | [RTS][RTS]               | [RTS][RTS], [Hyperdrive][Hyperdrive] |
| MWAX   | [Hyperdrive][Hyperdrive] | [Hyperdrive][Hyperdrive]             |

The links in the table will take you to the corresponding documentation.
Apart from [the MWA Telescope Wiki][Hyperdrive] (same link as given in the table), Hyperdrive also has some documentation on [its Github main page](https://github.com/MWATelescope/mwa_hyperdrive), and [its Github Wiki page](https://github.com/MWATelescope/mwa_hyperdrive/wiki).

The [RTS][RTS] link describes a workflow for preparing a calibration solution using the RTS.
This page describes the equivalent workflow for Hyperdrive solutions.
However, it should be noted that the visualisation tools used for Hyperdrive can also be used for RTS solutions, which are not documented there.

