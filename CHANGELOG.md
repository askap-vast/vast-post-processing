# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), with an added `List of PRs` section and links to the relevant PRs on the individal updates. This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased](https://github.com/askap-vast/vast-post-processing/compare/v2.1.1...HEAD)

### Added
- Added a filter that removes sources at the edges of observations when calculating the astromteric and flux corrections. These sources will be cropped out later. [#121](https://github.com/askap-vast/vast-post-processing/pull/121)
- Added a method that uses a linear fit with a Huber regressor to calculate the flux correction factor. [#123](https://github.com/askap-vast/vast-post-processing/pull/123)
- Added statsmodels to the poetry file.

### Changed
- Changed the way the flux offset factor is calculated, defaults to Huber regressor.

### Fixed
- Fixed error in extending source filtering that reversed numerator and denominator [#120](https://github.com/askap-vast/vast-post-processing/pull/120)
- Fixed how post-processing parameters were inserted in the catalogues/xml tables [#125](https://github.com/askap-vast/vast-post-processing/pull/125)

### Removed
- In crossmatch.py, removed the calculate_flux_offsets() function. Replaced with calculate_flux_offsets_median() [#122](https://github.com/askap-vast/vast-post-processing/pull/122)
- In crossmatch.py, removed the faulty median_abs_deviation() function. Replaced with astropy.stats.mad_std [#122](https://github.com/askap-vast/vast-post-processing/pull/122)
- [#122](https://github.com/askap-vast/vast-post-processing/pull/122): removed:  Changed the way the flux offset factor is calculated. Instead of using the slope of the linear fit to the flux_int vs. flux_int_reference plane, we now use the median of the distribution flux_int/flux_int_reference. Redundant now that we use the Huber regressor. Pull request was not implemented, closed.

### List of PRs
- [#125](https://github.com/askap-vast/vast-post-processing/pull/125): fixed: Fixed how post-processing parameters were inserted in the catalogues/xml tables.

- [#123](https://github.com/askap-vast/vast-post-processing/pull/123): added: Added a method that uses a linear fit with a Huber regressor to calculate the flux correction factor. Made this the default.

- [#122](https://github.com/askap-vast/vast-post-processing/pull/122): changed:  Changed the way the flux offset factor is calculated. Instead of using the slope of the linear fit to the flux_int vs. flux_int_reference plane, we now use the median of the distribution flux_int/flux_int_reference. This is similar to how the correction is calculate for the astrometry. Both flux and astrometry correction now use the median absolute deviation as an error. Pull request was not implemented, closed.

- [#121](https://github.com/askap-vast/vast-post-processing/pull/121): added: Added a filter that removes sources at the edges of observations when calculating the astromteric and flux corrections. These sources will be cropped out later.

- [#120](https://github.com/askap-vast/vast-post-processing/pull/120): fix: Fixed error in extending source filtering that reversed numerator and denominator

## [v2.1.1](https://github.com/askap-vast/vast-post-processing/releases/tag/v2.1.1) (2025-04-15)

### Fixed

- Fixed path image path matching for Stokes V processing [#117](https://github.com/askap-vast/vast-post-processing/pull/117)
- Fixed correction application for Stokes V processing [#117](https://github.com/askap-vast/vast-post-processing/pull/117)

### List of PRs

- [#117](https://github.com/askap-vast/vast-post-processing/pull/117): fix: Fixed Stokes V path finding and correction application

## [v2.1.0](https://github.com/askap-vast/vast-post-processing/releases/tag/v2.1.0) (2025-03-27)

### Added

- Allow users to specify whether or not to calculate Condon errors, and default to False [#114](https://github.com/askap-vast/vast-post-processing/pull/114)
- Added `reference_catalog` arg to Catalog class, which only sets the sbid, epoch and field for non-reference observations [#110](https://github.com/askap-vast/vast-post-processing/pull/110/)
- Allow users to specify quantisation level [#101](https://github.com/askap-vast/vast-post-processing/pull/101)
- Added changelog [#106](https://github.com/askap-vast/vast-post-processing/pull/106)

### Changed

- Changed compress_quant from 1024 to 16 [#113](https://github.com/askap-vast/vast-post-processing/pull/113)
- Changed crop size to 6.67 degrees [#111](https://github.com/askap-vast/vast-post-processing/pull/111)
- Change the flux scale corrections to use integrated rather than peak flux density [#109](https://github.com/askap-vast/vast-post-processing/pull/109)
- Generalised reference catalogue path [#108](https://github.com/askap-vast/vast-post-processing/pull/108)

### Fixed

- Fixed versioning issues [#115](https://github.com/askap-vast/vast-post-processing/pull/115)
- Fixed bug where critical error was being thrown if the reference images do not exist [#110](https://github.com/askap-vast/vast-post-processing/pull/110/)
- Fixed astrometric correction error - corrections were being applied in the wrong direction, resulting in the offset doubling [#104](https://github.com/askap-vast/vast-post-processing/pull/104)

### Removed

- Removed filtering on `has_siblings` as this should be already handled by the nearest-neighbour checks [#110](https://github.com/askap-vast/vast-post-processing/pull/110/)
- Removed hardcoded newest epoch variable [#107](https://github.com/askap-vast/vast-post-processing/pull/107)

### List of PRs

- [#115](https://github.com/askap-vast/vast-post-processing/pull/115): fix: Fixed versioning issues
- [#114](https://github.com/askap-vast/vast-post-processing/pull/114): feat: Allow users to specify whether or not to calculate Condon errors, and default to False
- [#113](https://github.com/askap-vast/vast-post-processing/pull/113): fix: Changed compress_quant from 1024 to 16
- [#111](https://github.com/askap-vast/vast-post-processing/pull/111): fix: Changed crop size to 6.67 degrees
- [#109](https://github.com/askap-vast/vast-post-processing/pull/109): feat: Switch flux scale corrections from peak to integrated
- [#110](https://github.com/askap-vast/vast-post-processing/pull/110/): feat, fix: add reference_catalog arg to catalog class, handle non-existence of reference images
- [#108](https://github.com/askap-vast/vast-post-processing/pull/108): feat: Generalised reference catalogue path
- [#104](https://github.com/askap-vast/vast-post-processing/pull/104): fix: Fix error in astrometric correction application
- [#101](https://github.com/askap-vast/vast-post-processing/pull/101): feat: Allow users to specify quantisation level
- [#106](https://github.com/askap-vast/vast-post-processing/pull/106): docs: Added changelog
- [#107](https://github.com/askap-vast/vast-post-processing/pull/107): fix: Removed hardcoded newest epoch variable

## [v2.0.0](https://github.com/askap-vast/vast-post-processing/releases/tag/v2.0.0) (2024-06-12)

This was the initial version 2 release that upgraded the code for use with the full VAST survey data. No changelog was kept.
