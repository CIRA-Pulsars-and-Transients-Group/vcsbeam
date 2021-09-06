# pfb\_filter

The files in this folder contain filter coefficients relating to various polyphase filter bank (PFB) operations.

## Description

The "raw" timeseries voltage data from the MWA pass through two distinct channelisation stages, before being delivered to the user via the Voltage Capture System (VCS; Tremblay et al. 2015). The first stage forms so-called "coarse channels", 1.28 MHz wide, at the Nyquist rate. The second stage forms the "fine channels", 10 kHz wide.

From late 2021 onwards, the legacy VCS is being replaced by the Phase 3 MWAX Correlator system, in which only the first-stage coarse channelisation is done.

The filters in this folder relate exclusively to the second stage channelisation.

## Contents

The coefficients for each filter described below can be found in this directory in files with the same name as the filter and the **.dat** extension.

- **FINEPFB**: The original integer coefficients used in the second-stage channlisation (the "fine PFB"). To be used as an *analysis* filter.
- **MIRROR**: A normalised (floating point) version of, but otherwise identical to, **FINEPFB.dat**. This is designed to be used as a *synthesis* filter
- **LSQ12**: An alternative *synthesis* filter, obtained via a least-squares optimisation method.

## More information

The above filters, and how they are determined, are described in more detail in the following paper:

McSweeney et al. (2020). *MWA tied-array processing III: Microsecond time resolution via a polyphase synthesis filter.* **PASA**, 37, e034. [https://doi.org/10.1017/pasa.2020.24](https://doi.org/10.1017/pasa.2020.24)
