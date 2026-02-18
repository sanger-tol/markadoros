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
markadoros build-database --preset bold --outdir db/ bold_release.fasta.gz

# 2. Search your reads
markadoros search-reads --index db/db.json --db BOLD_COI reads.fq.gz

# 3. Check results
less results.marker.result.tsv
```

## Detailed Usage

### Database Preparation

Use the `build-database` command to prepare marker gene sequences for searching:

```bash
markadoros build-database --preset bold --outdir db/ /path/to/bold/release.fasta.gz
```

**Available presets:**

- `bold` - [BOLD Systems](https://www.boldsystems.org/) releases
- `unite` - [UNITE](https://unite.ut.ee/) ITS database

**Custom parameters:**

```bash
markadoros build-database --db-params params.json --outdir db/ sequences.fasta
```

The output includes:
- `db.json` - Index of available databases and their parameters
- `<database_name>/db` - MMSeqs2 database files

### Searching Reads

Use the `search-reads` command to identify barcode genes:

```bash
markadoros search-reads --index db/db.json reads.fq.gz
```

**Options:**
- `--db <name>` - Search specific database (default: all)
- `--nreads <N>` - Limit to first N reads
- `--prefix <name>` - Output file prefix (default: input filename)
- `--platform illumina|illumina_rnaseq|pacbio_hifi|pb|oxford_nanopore|ont` - sequencing platform of input reads

**Output files:**
- `<prefix>.marker.contigs.fa` - Assembled marker gene sequences
- `<prefix>.marker.result.tsv` - Search results (BLAST tabular format with extensions)

### Database Parameters JSON

Create a custom `params.json` for your marker genes:

```json
{
    "parameters": {
        "deduplicate": true,
        "seq_id_regex": "^([^|]+)",
        "marker_id_regex": "^[^|]+\\|([^|]+)\\|",
        "taxon_id_regex": "\\|[^|]*?([^,None][^,]*)(?:,None)*$"
    },
    "databases": {
        "MY_MARKER": {
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
- `seq_id_regex` - Extract sequence identifier from header (must have capture group)
- `marker_id_regex` - Extract marker name (can be null for single-marker databases)
- `taxon_id_regex` - Extract taxonomic classification

**Database options:**
- `marker` - Marker gene name
- `min_length` - Minimum sequence length to include
- `min_seq_id` - Minimum sequence identity for search hits
- `min_aln_len` - Minimum alignment length for search hits

## Author

- Jim Downie
