# hifi_trimmer: changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.2.1](https://github.com/sanger-tol/markadoros/releases/tag/v1.2.1)] - [2026-08-20]

### Fixed

- Update dependencies and remove extreaneous scikit-learn dependency, which was previously required for pymmseqs.

## [[1.2.0](https://github.com/sanger-tol/markadoros/releases/tag/v1.2.0)] - [2026-08-10]

### Changed

- The results are no longer sorted by coverage, so the top result is the top result across all contigs rather than the contig with the highest coverage.
- When building a COI database from a BOLD TSV file, the taxonomic makeup of each BIN is now stored as part of the database. This is then attached to each result so that other taxon names within the bin can be identified in cases of multi-taxon BINs.
- The `is_expected_taxon` field in the summary JSON is now a simple boolean value. An additional field, `expected_taxon_match_type`, now provides the match type for the expected taxon (one of `exact`, `synonym`, `bold_bin_contains_expected_taxon`, `bold_bin_contains_synonym`).

## Fixed

- When using `--save-contigs`, if the pipeline has fallen back to searching reads these will now be correctly saved with a `.gz` extension.

## [[1.1.0](https://github.com/sanger-tol/markadoros/releases/tag/v1.1.0)] - [2026-04-23]

### Changed

- Markadoros now reports the top hit and the top hit for the expected taxon separately in the summary JSON. 

## [[1.0.0](https://github.com/sanger-tol/markadoros/releases/tag/v1.0.0)] - [2026-04-23]

### Added

- Added options `--find-goat-synonyms` to `markadoros search` to get a list of synonyms from the expected taxon from [GoaT](https://goat.genomehubs.org), 
  `--synonyms` to provide a comma-separated list of synonyms.
- `markadoros search` will now fall back and search the first `--fallback-reads` reads assembly fails but reads containing markers were found.

### Fixed

- If an expected taxon is provided, the results are now correctly sorted by bit score rather than fident and alnlen.

## [[0.1.0](https://github.com/sanger-tol/markadoros/releases/tag/v0.100)] - [2026-04-01]

Initial release of markadoros.
