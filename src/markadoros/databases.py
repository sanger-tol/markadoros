import gzip
import hashlib
import io
import json
import re
import shutil
from contextlib import redirect_stdout
from pathlib import Path

import click
import pysam
from Bio.Seq import Seq
from jsonschema import ValidationError, validate
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig


class DatabaseCreator:
    def __init__(
        self,
        outdir: Path,
        db_dict: dict,
    ) -> None:
        ## Stage output dir
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)

        ## Tempdir staging
        self.rmtmp = False
        self.tmpdir = self.outdir / "tmp"
        self.tmpdir.mkdir(parents=True, exist_ok=True)

        ## Setup split dict
        self.split_dict = db_dict

        ## Extract database params
        self.deduplicate = db_dict.get("parameters").get("deduplicate")
        self.seq_id_regex = db_dict.get("parameters").get("seq_id_regex")
        self.marker_id_regex = db_dict.get("parameters").get("marker_id_regex")
        self.taxon_id_regex = db_dict.get("parameters").get("taxon_id_regex")

        self.databases = db_dict.get("databases")

    def _get_canonical_sequence(self, seq: str) -> str:
        """
        Get the canonical form of a sequence (lexicographically smaller of seq and its reverse complement).
        """
        seq_str = str(seq).upper()
        rev_comp = str(Seq(seq_str).reverse_complement())

        return min(seq_str, rev_comp)

    def _compute_sequence_hash(self, sequence: str) -> str:
        """Compute MD5 hash of canonical sequence."""
        canonical_seq = self._get_canonical_sequence(sequence)
        return hashlib.md5(canonical_seq.encode()).hexdigest()

    def _write_fasta_record(self, handle, name, seq) -> None:
        """Write a FASTA record as gzip-compressed bytes."""
        handle.write(b">")
        handle.write(name.encode() + b"\n")
        handle.write(seq.encode() + b"\n")

    def _process_fasta_file(
        self,
        fasta: Path,
        output_handles: dict,
        seen_sequences: dict = None,
    ) -> None:
        """
        Generic FASTA file processor with progress tracking and deduplication.

        Args:
            fasta: Path to input FASTA file
            output_handles: Dict mapping db_name -> gzip file handle
            seen_sequences: Dict mapping db_name -> set() for dedup tracking
        """
        with pysam.FastxFile(fasta, persist=False) as fh:
            for record in fh:
                seq_hash = None
                if self.deduplicate:
                    seq_hash = self._compute_sequence_hash(record.sequence)

                yield record, seq_hash

    def _process_fasta(self, fasta: Path) -> dict:
        """Process a FASTA file, splitting by markers with optional deduplication."""

        ## Extract the list of marker genes to check against.
        markers = [info.get("marker") for info in self.databases.values()]

        splitting = True
        marker = None
        if len(markers) == 1:
            splitting = False
            marker = markers[0]

        click.echo(f"Splitting {fasta.name} by markers: {', '.join(markers)}")

        ## Open output files for each database
        output_handles = {
            db: gzip.open(self.tmpdir / f"{db}.fa.gz", "ab")
            for db in self.databases.keys()
        }

        ## Track which sequences are seen if we are deduplicatng
        seen_sequences = (
            {db: set() for db in self.databases.keys()} if self.deduplicate else None
        )

        ## Count records processed per database
        record_counts = {db: 0 for db in self.databases.keys()}

        try:
            seq_count = 0
            for record, seq_hash in self._process_fasta_file(
                fasta, output_handles, seen_sequences
            ):
                seq_count += 1
                if seq_count % 100000 == 0:
                    click.echo(f"Processed {seq_count:,} sequences", err=True)

                if record.comment is None:
                    header = record.name
                else:
                    header = "_".join([record.name, record.comment])

                seq_id = re.search(self.seq_id_regex, header).group(1)
                taxon_id = re.search(self.taxon_id_regex, header).group(1)

                if not splitting and self.marker_id_regex is None:
                    marker_id = marker
                else:
                    marker_id = re.search(self.marker_id_regex, header).group(1)

                output_name = "|".join([seq_id, marker_id, taxon_id])
                if len(output_name) == 0:
                    raise ValueError(
                        f"Could not find any matches in FASTA header: {header}"
                    )

                for db_name, marker_info in self.databases.items():
                    if splitting and marker_info.get("marker") not in marker_id:
                        continue

                    if len(record.sequence) < marker_info.get("min_length"):
                        continue

                    if self.deduplicate and seq_hash in seen_sequences[db_name]:
                        continue

                    if self.deduplicate:
                        seen_sequences[db_name].add(seq_hash)

                    record_counts[db_name] += 1
                    self._write_fasta_record(
                        output_handles[db_name], output_name, record.sequence
                    )

        except IOError as e:
            raise IOError(f"Error processing {fasta.name}: {e}") from e

        finally:
            for handle in output_handles.values():
                handle.close()

        for db, count in record_counts.items():
            if count == 0:
                click.echo(f"No records found for {db}")

        return {
            db: {**self.databases[db], "fasta": (self.tmpdir / f"{db}.fa.gz").resolve()}
            for db in self.databases.keys()
            if record_counts[db] > 0
        }

    def _append_to_db_index(self, marker_name: str, db_entry: dict) -> None:
        """Append marker and database information to JSON index file.

        Args:
            marker_name: Name of the marker/database
            db_entry: Dictionary containing all params and paths (fasta, db, marker, min_length, etc.)
        """
        db_index_path = self.outdir / "db.json"

        try:
            with open(db_index_path, "r") as f:
                db_dict = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            db_dict = {}

        db_dict[marker_name] = {
            **db_entry,
        }

        with open(db_index_path, "w") as f:
            json.dump(db_dict, f, indent=2)

    def create_marker_database(
        self,
        fasta: Path,
    ):
        """
        Create MMSeqs2 databases. If is_bold_fasta is true, we assume we have a
        BOLD FASTA file and split it by markers prior to DB creation.
        """
        outfiles_dict = self._process_fasta(fasta=fasta)

        for database, params in outfiles_dict.items():
            click.echo(f"Building MMSeqs2 database for {database}... ", nl=False)

            db_path = (self.outdir / database / "db").resolve()

            with redirect_stdout(io.StringIO()):
                db_config = CreateMMSeqsDBConfig(
                    fasta_file=params.get("fasta"), sequence_db=db_path, dbtype=2, v=1
                )
                db_config.run()

            # Add the MMSeqs DB path to params before storing
            db_entry = {**params, "db": str(db_path)}
            db_entry.pop("fasta", None)
            self._append_to_db_index(database, db_entry)

            click.echo("done")

        if self.rmtmp:
            shutil.rmtree(self.tmpdir)


def validate_db(database: Path) -> dict:
    """
    Validate a database JSON file against the expected schema.

    Args:
        database: Path to the database.json file

    Returns:
        dict: The validated database configuration

    Raises:
        ValueError: If the file doesn't exist, is invalid JSON, or doesn't conform to schema
        ValidationError: If the JSON doesn't match the expected schema
    """
    if not database.exists():
        raise ValueError(f"{database} is not a valid database!")

    try:
        with open(database, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {database}: {e}")

    db_schema = {
        "type": "object",
        "additionalProperties": {
            "type": "object",
            "properties": {
                "marker": {"type": "string"},
                "min_length": {"type": "integer", "minimum": 1},
                "min_seq_id": {"type": "number", "minimum": 0, "maximum": 1},
                "min_aln_len": {"type": "integer", "minimum": 1},
                "db": {"type": "string"},
            },
            "required": ["marker", "min_length", "min_seq_id", "min_aln_len", "db"],
        },
    }

    try:
        validate(instance=data, schema=db_schema)
    except ValidationError as e:
        raise ValueError(
            f"Database configuration in {database} does not match expected schema: {e.message}"
        ) from e

    return data
