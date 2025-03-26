# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), with an added `List of PRs` section and links to the relevant PRs on the individal updates. This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased](https://github.com/askap-vast/vast-post-processing/compare/v2.0.0...HEAD)

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

- Fixed bug where critical error was being thrown if the reference images do not exist [#110](https://github.com/askap-vast/vast-post-processing/pull/110/)
- Fixed astrometric correction error - corrections were being applied in the wrong direction, resulting in the offset doubling [#104](https://github.com/askap-vast/vast-post-processing/pull/104)

### Removed

- Removed filtering on `has_siblings` as this should be already handled by the nearest-neighbour checks [#110](https://github.com/askap-vast/vast-post-processing/pull/110/)
- Removed hardcoded newest epoch variable [#107](https://github.com/askap-vast/vast-post-processing/pull/107)

### List of PRs

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
