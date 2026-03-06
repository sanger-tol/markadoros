# Markadoros - fast barcode assembly and identification from raw sequencing data

## Introduction

Markadoros is a Python tool for the identification and assembly of barcode genes in raw sequencing data, using
[MMSeqs2](https://github.com/soedinglab/MMseqs2) for searching and either [Spades](https://github.com/ablab/spades) or [Hifiasm](https://github.com/chhylp123/hifiasm) for assembly. Using an MMSeqs2 database of marker gene sequences, markadoros searches a set of input reads (providable in either FASTX or CRAM format) to quickly pre-filter reads that map to the target marker gene. Reads that match any barcode sequence in the database are extracted and assembled using an appropriate assembler. The resulting assembled contigs are then searched again using the same database, using more sensitive thresholds, to identify the marker genes from the input dataset.

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
# 1. Prepare a database from BOLD release
markadoros database --preset bold --outdir db/ bold_release.fasta.gz

# 2. Search your reads
markadoros search --type illumina --index db/db.json --db BOLD_COI reads.fq.gz

# 3. Check results
less reads.BOLD_COI.summary.json
```

## Detailed Usage

### Database preparation

Use the `database` command to prepare marker gene sequences for searching:

```bash
markadoros database --preset bold --outdir db/ /path/to/bold/release.fasta.gz
```

You can choose a preset to process some pre-defined input FASTA releases correctly. Currently, there are two supported presets:

- `bold` - [BOLD Systems](https://www.boldsystems.org/) [general FASTA release](https://bench.boldsystems.org/index.php/datapackages/Latest)
- `unite` - [UNITE](https://unite.ut.ee/) [general FASTA release](https://unite.ut.ee/repository.php)

If your FASTA release does not conform to the above, you can build a custom database. You will
need to create a params JSON file for your database (see below), and the input FASTA headers must be of the format:

```
><unique_id>|<marker>|<taxon_name>|<lineage or null>
```

You can then process the input as follows:

```bash
markadoros build-database --db-params params.json --outdir db/ sequences.fasta
```

The output of DB preparation includes:

- `db.json` - Index of available databases and their parameters
- `<database_name>/db*` - MMSeqs2 database files

If you build a second database pointing to the same output directory, the existing index file will be updated to include the new index.

### Searching for barcodes 

Use the `search` command to identify barcode genes in a set of reads or a set of contigs:

```bash
markadoros search --type illumina --index db/db.json reads.fq.gz
```

**Options:**
- `--type, -t` - Input data type (required): `[sr|short|illumina]`, `rnaseq`, `[pacbio_hifi|pb]`, `[oxford_nanopore|ont]`, or `contigs`
- `--index, -i` - Path to database index JSON (required)
- `--db <name>` - Search specific database (default: all)
- `--nreads, -n <N>` - Limit to first N reads
- `--prefix, -p <name>` - Output file prefix (default: input filename)
- `--outdir, -o <path>` - Output directory (default: current directory)
- `--threads, -t <N>` - Number of threads (default: 1)
- `--expected_taxon <name>` - Expected taxon binomial name for validation
- `--cleanup/--no-cleanup` - Clean up temporary files (default: cleanup)

**Output files:**
- `<prefix>.marker.contigs.fa` - Assembled marker gene sequences
- `<prefix>.marker.result.json` - Results in JSON format

### Database Parameters JSON

Create a custom `params.json` for your marker genes:

```json
{
    "parameters": {
        "deduplicate": true
    },
    "databases": {
        "MY_MARKER_DATABASE": {
            "marker": "COI-5P",
            "min_length": 200,
            "min_seq_id": 0.96,
            "min_aln_len": 450
        }
    }
}
```

**Parameter descriptions:**
- `deduplicate` - Remove duplicate sequences

**Database options:**
- `marker` - Marker gene name
- `min_length` - Minimum sequence length to include
- `min_seq_id` - Minimum sequence identity for search hits
- `min_aln_len` - Minimum alignment length for search hits

## Author

- Jim Downie
