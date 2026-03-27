# Markadoros - fast barcode assembly and identification from raw sequencing data

## Introduction

Markadoros is a Python tool for the identification and assembly of barcode genes in raw sequencing data, using
[MMSeqs2](https://github.com/soedinglab/MMseqs2) for searching and either [SPAdes](https://github.com/ablab/spades) or [Hifiasm](https://github.com/chhylp123/hifiasm) for assembly. Using an MMSeqs2 database of marker gene sequences, markadoros searches a set of input reads (providable in either FASTX or CRAM format) to quickly pre-filter reads that map to the target marker gene. Reads that match any barcode sequence in the database are extracted and assembled using an appropriate assembler. The resulting assembled contigs are then searched again using the same database, using more sensitive thresholds, to identify the marker genes from the input dataset.

The name *markadoros* comes from the Greek word for a marker pen, *μαρκαδόρος*.

## Requirements

### System Dependencies

You will need the following external tools installed to run Markadoros:

- [MMSeqs2](https://github.com/soedinglab/MMseqs2)
- [SPAdes](https://github.com/ablab/spades)
- [Hifiasm](https://github.com/chhylp123/hifiasm)

### Python Dependencies

Markadoros requires Python 3.11+ and the following packages:

- biopython (≥1.86)
- click (≥8.3.1)
- jsonschema (≥4.26.0)
- pandas (≥2.3.3)
- pymmseqs (≥1.0.5)
- pysam (≥0.23.3)
- scikit-learn (≥1.8.0)

These are automatically installed when you install Markadoros.

## Installation

### From source

```bash
pip install -e .
```

### Dependencies

Install the required external tools with Conda:

```bash
conda install -c bioconda mmseqs2 spades hifiasm
```

## Quick Start

```bash
# 1. Build a database from a BOLD release
markadoros database -x bold --marker COI --marker rbcL --prefix BOLD --outdir db/ <bold_release.fasta.gz>

# 2. Search your reads
markadoros search -x illumina --index db/db.json --db BOLD.COI <reads.fq.gz>

# 3. Check results
less reads.COI.summary.json
```

## Detailed Usage

### Database Preparation

Use the `database` command to prepare marker gene sequences for searching:

```bash
markadoros database -x bold \
    --marker COI --marker rbcL --marker CYTB \
    --marker matK --marker 18S --marker 28S \
    --prefix BOLD \
    --outdir db/ \
    /path/to/bold/release.fasta.gz
```

**Required options:**

- `--marker <name>` - Marker gene name to extract from input FASTA. Can be specified multiple times for FASTAs containing multiple markers (e.g., `--marker COI-5P --marker ITS`). Inexact matches are allowed. Required unless `--header-type` is `unite`.
- `--prefix <name>` - Prefix for output database names

**Additional options:**

- `-x, --header-type <type>` - Use a preset header processor: `bold`, `unite`, `silva_lsu`, `silva_ssu`.
- `--min-length <N>` - Minimum sequence length to retain (default: 200)
- `--deduplicate/--no-deduplicate` - Deduplicate identical sequences. The first record for each identical sequence is kept. (default: deduplicate)
- `--cluster/--no-cluster` - Cluster sequences using MMSeqs2's linear clustering algorithm (default: no cluster)
- `--cluster_min_seq_id` - Cluster sequences at this percentage identity threshold (default: 0.99)
- `--cluster_coverage` - Overlap between two sequences required for clustering (default: 0.8)
- `--create-index/--no-create-index` - Create an MMSeqs2 index for each marker database. This may improve speed for larger databases, but can cause IO issues if multiple processes access the same database. (default: False)
- `--exclude-file` - New line-separated list of regular expressions. If a header matches a regular expression it will be skipped.
- `-o, --outdir <path>` - Output directory (default: `./markadoros.db`)
- `--cleanup/--no-cleanup` - Clean up temporary files (default: cleanup)
- `-t, --threads <N>` - Number of threads for MMSeqs2 (default: 1)

**Header types:**

- `bold` - [BOLD Systems](https://www.boldsystems.org/) [general FASTA release](https://bench.boldsystems.org/index.php/datapackages/Latest)
- `unite` - [UNITE](https://unite.ut.ee/) [general FASTA release](https://unite.ut.ee/repository.php)
- `silva_lsu` - [SILVA](https://www.arb-silva.de/) [LSU Release](https://www.arb-silva.de/current-release/Exports)
- `silva_ssu` - [SILVA](https://www.arb-silva.de/) [SSU Release](https://www.arb-silva.de/current-release/Exports)
- 
If your FASTA release does not conform to the above presets, omit the `-x` option. Your input FASTA headers must then be formatted as:

```
><unique_id>|<marker>|<taxon_name>|<lineage>
```

**Output files:**

- `db.json` - Index of available databases and their parameters
- `<prefix>_<marker>/db*` - MMSeqs2 database files
- `<prefix>_<marker>/taxon.json.gz` - JSON file counting the number of available sequences per taxon.

If you build additional databases pointing to the same output directory, the existing index file will be updated to include the new entries.

### Searching for Barcodes

Use the `search` command to identify barcode genes in a set of reads or pre-assembled contigs:

```bash
markadoros search -x illumina --index db/db.json reads.fq.gz
```

**Required options:**

- `-x, --type <type>` - Input data type (see table below)
- `-i, --index <path>` - Path to database index JSON

**Input type aliases:**

| Platform               | Accepted values                      |
|------------------------|--------------------------------------|
| Illumina / short reads | `sr`, `short`, `illumina`            |
| Short read RNA-seq     | `rnaseq`                             |
| PacBio HiFi            | `pb`, `pacbio`, `pacbio_hifi`        |
| Oxford Nanopore        | `ont`, `nanopore`, `oxford_nanopore` |
| Pre-assembled contigs  | `contigs`                            |

**Additional options:**

- `--db <name>` - Search a specific database only (default: search all databases in index)
- `--expected_taxon <name>` - Expected taxon binomial name for validation
- `-n, --nreads <N>` - Limit to first N reads
- `-m, --min_seq_id <float>` - Minimum sequence identity for hits (default: 0.96)
- `-l, --min_aln_len <int>` - Minimum alignment length for hits (default: 450)
- `--cleanup/--no-cleanup` - Clean up temporary files (default: cleanup)
- `--db-to-tmpdir/--no-db-to-tmpdir` - Temporarily copy the database to the temporary directory. Can reduce IO and improve speed if multiple processes would access the database simultaneously. Not suggested if the databases include indexes. (default: True)
- `-t, --threads <N>` - Number of threads (default: 1)
- `-p, --prefix <name>` - Output file prefix (default: input filename)
- `-o, --outdir <path>` - Output directory (default: current directory)

**Output files:**

- `<prefix>.<marker>.summary.json` - Results in JSON format

### Output JSON file

The search subtool outputs a JSON file summarising the search. It has the following format:

```json
{
  "input": {
      "file": <path>, // input file path
      "n_reads": <int>, // number of reads searched
      "n_aligned_reads": <int>, // number of reads aligned to database
      "marker": <string>, // marker gene in database
      "database": <path>, // path to mmseqs database
      "contig_stats": <dict>, // number, total length and n50 of assembled contigs
      "expected_taxon": <dict>, // name of asserted taxon and number of records present in database
  },
  "summary": {
      "n_contigs_with_hits": <int>, // number of assembled contigs with search hits
      "n_expected_taxon_hits": <int>, // number of hits for the asserted taxon
      "top_result": <dict>, // top result (highest bitscore) for asserted taxon if found or overall if not
      "taxon_summary": <dict>, // per-taxon summary of results - number of hits, min and max %ID and alignment lengths, and sequence of top hit
  },
  "results": <list>, // for each contig with results, a list of all hits
  "run_info": <dict> // tool information
}

```

### Example Workflows

**Building a marker database for the BOLD FASTA release:**

```bash
markadoros database \
    --header-type bold \
    --outdir db/ \
    --marker COI --marker rbcL --marker CYTB \
    --marker matK --marker 18S --marker 28S \
    --prefix BOLD \
    --cluster \
    --no-deduplicate \
    --threads 16 \
    BOLD_Public.20-Feb-2026.fasta.gz
```

**Searching PacBio HiFi reads with expected taxon:**

```bash
markadoros search -x pb \
    --index db/db.json \
    --db BOLD.COI \
    --expected_taxon "Halyzia sedecimguttata" \
    --threads 16 \
    --nreads 20000 \
    pacbio_reads.fasta.gz
```

## Author

- Jim Downie
