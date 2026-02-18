import gzip
import hashlib
import re
from pathlib import Path

import click
import pysam
from Bio.Seq import Seq


class FASTAProcessor:
    """Handles FASTA file reading, header parsing, and deduplication."""

    def __init__(
        self,
        seq_id_regex: str,
        marker_id_regex: str,
        taxon_id_regex: str,
        databases: dict,
        deduplicate: bool = False,
        tmpdir: Path = None,
    ):
        self.seq_id_regex = re.compile(seq_id_regex)
        self.marker_id_regex = re.compile(marker_id_regex) if marker_id_regex else None
        self.taxon_id_regex = re.compile(taxon_id_regex)
        self._databases = databases
        self._deduplicate = deduplicate
        self._tmpdir = tmpdir or Path.cwd()

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

    def _extract_header_fields(self, record) -> tuple:
        """Extract seq_id, marker_id, and taxon_id from a FASTA record header."""
        if record.comment is None:
            header = record.name
        else:
            header = "_".join([record.name, record.comment])

        try:
            seq_id = self.seq_id_regex.search(header).group(1)
            taxon_id = self.taxon_id_regex.search(header).group(1)

            # Determine marker_id based on whether we're splitting
            markers = [info.get("marker") for info in self._databases.values()]
            if len(markers) == 1 and self.marker_id_regex is None:
                marker_id = markers[0]
            else:
                marker_id = self.marker_id_regex.search(header).group(1)

            return seq_id, marker_id, taxon_id

        except AttributeError as e:
            raise ValueError(
                f"Could not extract fields from FASTA header: {header}"
            ) from e

    def _should_include_record(
        self,
        record,
        marker_id: str,
        seq_hash: str,
        db_name: str,
        marker_info: dict,
        seen_sequences: dict,
    ) -> bool:
        """Check if a record should be included in the database."""
        # Check if marker matches
        markers = [info.get("marker") for info in self._databases.values()]
        if len(markers) > 1 and marker_info.get("marker") not in marker_id:
            return False

        # Check minimum length
        if len(record.sequence) < marker_info.get("min_length"):
            return False

        # Check for duplicates
        if self._deduplicate and seq_hash in seen_sequences[db_name]:
            return False

        return True

    def _write_fasta_record(self, handle, name: str, seq: str) -> None:
        """Write a FASTA record to a gzip file handle."""
        handle.write(b">")
        handle.write(name.encode() + b"\n")
        handle.write(seq.encode() + b"\n")

    def process(self, fasta: Path) -> dict:
        """Process a FASTA file, splitting by markers with optional deduplication."""
        markers = [info.get("marker") for info in self._databases.values()]
        click.echo(f"Splitting {fasta.name} by markers: {', '.join(markers)}")

        # Open output files for each database
        output_handles = {
            db: gzip.open(self._tmpdir / f"{db}.fa.gz", "ab")
            for db in self._databases.keys()
        }

        # Track which sequences are seen if we are deduplicating
        seen_sequences = (
            {db: set() for db in self._databases.keys()} if self._deduplicate else None
        )

        # Count records processed per database
        record_counts = {db: 0 for db in self._databases.keys()}

        try:
            seq_count = 0
            with pysam.FastxFile(fasta, persist=False) as fh:
                for record in fh:
                    seq_count += 1
                    if seq_count % 100000 == 0:
                        click.echo(f"Processed {seq_count:,} sequences", err=True)

                    # Extract header fields
                    seq_id, marker_id, taxon_id = self._extract_header_fields(record)
                    output_name = "|".join([seq_id, marker_id, taxon_id])

                    # Compute sequence hash for deduplication
                    seq_hash = None
                    if self._deduplicate:
                        seq_hash = self._compute_sequence_hash(record.sequence)

                    # Write to appropriate database files
                    for db_name, marker_info in self._databases.items():
                        if not self._should_include_record(
                            record,
                            marker_id,
                            seq_hash,
                            db_name,
                            marker_info,
                            seen_sequences,
                        ):
                            continue

                        if self._deduplicate:
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

        # Log databases with no records
        for db, count in record_counts.items():
            if count == 0:
                click.echo(f"No records found for {db}")

        # Build output dictionary with only databases that have records
        return {
            db: {
                **self._databases[db],
                "processed_fasta": (self._tmpdir / f"{db}.fa.gz").resolve(),
                "n_seqs": record_counts[db],
            }
            for db in self._databases.keys()
            if record_counts[db] > 0
        }
