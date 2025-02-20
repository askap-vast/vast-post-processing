# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), with an added `List of PRs` section and links to the relevant PRs on the individal updates. This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

##[Unreleased](https://github.com/askap-vast/vast-post-processing/compare/v2.0.0...HEAD)

### Added

- Allow users to specify quantisation level [#101](https://github.com/askap-vast/vast-post-processing/pull/101)
- Added changelog [#106](https://github.com/askap-vast/vast-post-processing/pull/106)

### Changed

### Fixed

- Fixed astrometric correction error - corrections were being applied in the wrong direction, resulting in the offset doubling [#104](https://github.com/askap-vast/vast-post-processing/pull/104)

### Removed

- Removed hardcoded newest epoch variable [#107](https://github.com/askap-vast/vast-post-processing/pull/107)

### List of PRs

- [#104](https://github.com/askap-vast/vast-post-processing/pull/104): fix: Fix error in astrometric correction application
- [#101](https://github.com/askap-vast/vast-post-processing/pull/101): feat: Allow users to specify quantisation level
- [#106](https://github.com/askap-vast/vast-post-processing/pull/106): docs: Added changelog
- [#107](https://github.com/askap-vast/vast-post-processing/pull/107): fix: Removed hardcoded newest epoch variable

## [v2.0.0](https://github.com/askap-vast/vast-post-processing/releases/tag/v2.0.0) (2024-06-12)

This was the initial version 2 release that upgraded the code for use with the full VAST survey data. No changelog was kept.
