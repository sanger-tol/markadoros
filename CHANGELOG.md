# hifi_trimmer: changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0](https://github.com/sanger-tol/markadoros/releases/tag/v1.0.0)] - [2026-04-23]

### Added

- Added options `--find-goat-synonyms` to `markadoros search` to get a list of synonyms from the expected taxon from [GoaT](https://goat.genomehubs.org), 
  `--synonyms` to provide a comma-separated list of synonyms.
- `markadoros search` will now fall back and search the first `--fallback-reads` reads assembly fails but reads containing markers were found.

### Fixed

- If an expected taxon is provided, the results are now correctly sorted by bit score rather than fident and alnlen.

## [0.1.0](https://github.com/sanger-tol/markadoros/releases/tag/v0.100)] - [2026-04-01]

Initial release of markadoros.
