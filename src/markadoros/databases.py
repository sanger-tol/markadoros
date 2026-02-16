import gzip
import hashlib
import json
import shutil
from pathlib import Path

import click
import pysam
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig

from markadoros.constants import VALID_MARKERS
from markadoros.utils import get_canonical_sequence


class DatabaseCreator:
    def __init__(self, outdir: Path, min_length: int, deduplicate: bool) -> None:
        self.outdir = outdir
        self.min_length = min_length
        self.deduplicate = deduplicate

    def _split_barcodes(
        self,
        fasta: Path,
        marker_dict: dict[str, str],
        tmpdir: Path,
    ) -> dict:
        """Split FASTA records by marker type, with optional deduplication."""
        click.echo(
            f"Splitting {fasta.name} by markers: {', '.join(marker_dict.keys())}"
        )

        # Initialize output handles and deduplication tracking
        output_handles = {
            marker: gzip.open(tmpdir / f"{marker}.fa.gz", "ab")
            for marker in marker_dict.keys()
        }
        seen_sequences = (
            {marker: set() for marker in marker_dict.keys()}
            if self.deduplicate
            else None
        )

        # Build a lookup dict for faster marker matching
        marker_search_map = {
            search_string: marker_name
            for marker_name, search_string in marker_dict.items()
        }

        try:
            seq_count = 0
            with pysam.FastxFile(fasta, persist=False) as fh:
                for record in fh:
                    seq_count += 1
                    if seq_count % 100000 == 0:
                        click.echo(f"Processed {seq_count:,} sequences", err=True)

                    seq_len = len(record.sequence)
                    if seq_len < self.min_length:
                        continue

                    # Parse marker once
                    try:
                        record_marker = record.name.split("|")[1]
                    except IndexError:
                        continue

                    # Compute hash once if deduplicating
                    seq_hash = None
                    if self.deduplicate:
                        canonical_seq = get_canonical_sequence(record.sequence)
                        seq_hash = hashlib.md5(canonical_seq.encode()).hexdigest()

                    # Find matching marker
                    for search_string, marker_name in marker_search_map.items():
                        if search_string not in record_marker:
                            continue

                        if self.deduplicate and seq_hash in seen_sequences[marker_name]:
                            continue

                        if self.deduplicate:
                            seen_sequences[marker_name].add(seq_hash)

                        # Write directly as bytes (avoid intermediate string formatting)
                        handle = output_handles[marker_name]
                        handle.write(b">")
                        handle.write(record.name.encode() + b"\n")
                        handle.write(record.sequence.encode() + b"\n")

        except IOError as e:
            raise IOError(f"Could not open file: {e}") from e

        finally:
            for handle in output_handles.values():
                handle.close()

        return {marker: tmpdir / f"{marker}.fa.gz" for marker in marker_dict.keys()}

    def _append_to_db_index(self, marker_name: str, db_path: str) -> None:
        """Append marker and database path to JSON index file."""
        db_index_path = self.outdir / "db.json"

        try:
            with open(db_index_path, "r") as f:
                db_dict = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            db_dict = {}

        db_dict[marker_name] = db_path

        with open(db_index_path, "w") as f:
            json.dump(db_dict, f, indent=2)

    def create_marker_database(
        self, fasta: Path, marker: str | None, is_bold_fasta: bool = False
    ):
        """
        Create MMSeqs2 databases. If is_bold_fasta is true, we assume we have a
        BOLD FASTA file and split it by markers prior to DB creation.
        """
        if is_bold_fasta:
            tmpdir = self.outdir / "tmp"
            tmpdir.mkdir(parents=True, exist_ok=True)

            outfiles_dict = self._split_barcodes(
                fasta=fasta, marker_dict=VALID_MARKERS, tmpdir=tmpdir
            )

            for marker_name, fasta_file in outfiles_dict.items():
                click.echo(f"Building MMSeqs2 database for {marker_name}")

                db_path = f"{self.outdir}/{marker_name}/db"

                db_config = CreateMMSeqsDBConfig(
                    fasta_file=fasta_file, sequence_db=db_path, db_type=2, v=1
                )
                db_config.run()

                self._append_to_db_index(marker_name, db_path)

            shutil.rmtree(tmpdir)
        else:
            click.echo(f"Building {marker} MMSeqs2 database from {fasta.name}")

            db_path = self.outdir / marker / "db"

            db_config = CreateMMSeqsDBConfig(
                fasta_file=fasta, sequence_db=db_path, db_type=2, v=1
            )
            db_config.run()

            self._append_to_db_index(marker, db_path)


def validate_db(database: Path) -> Path:
    """
    Validate a database to check it is, or has, a JSON file with correct schema.

    Expected schema: {"marker_name": "path/to/db", ...}
    """
    if database.is_dir():
        json_file = database / "db.json"
    else:
        json_file = database

    if not json_file.exists():
        raise ValueError(f"{database} is not a valid database!")

    try:
        with open(json_file, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {json_file}: {e}")

    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object at root level in {json_file}")

    for key, value in data.items():
        if not isinstance(key, str):
            raise ValueError(
                f"Expected marker names (keys) to be strings in {json_file}"
            )
        if not isinstance(value, str):
            raise ValueError(
                f"Expected database paths (values) to be strings in {json_file}"
            )

    return data
