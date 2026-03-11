import hashlib
import os
from pathlib import Path
from typing import Callable

import pgzip
import pysam
from Bio.Seq import Seq
from loguru import logger

from markadoros.header_processor import process_generic_header


class DatabaseFASTAProcessor:
    """Handles FASTA file reading, header parsing, and deduplication."""

    def __init__(
        self,
        deduplicate: bool = True,
        header_processor: Callable[[str], tuple[str, str, str]] | None = None,
        tmpdir: Path | None = None,
        min_length: int = 200,
        threads: int = 0,
    ):
        self._deduplicate = deduplicate
        self._min_length = min_length
        self._tmpdir = tmpdir or Path.cwd()
        # Use all available CPUs if threads is 0
        self._threads = threads if threads > 0 else os.cpu_count() or 4

        if header_processor is None:
            self.header_processor = process_generic_header
        else:
            self.header_processor = header_processor

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

    def _marker_matches(self, marker: str, db_marker: str) -> bool:
        """Check if a marker matches the database marker."""
        return marker.startswith(db_marker)

    def _is_duplicate(self, seq_hash: str, db_seen_sequences: set[str] | None) -> bool:
        """Check if a sequence hash has already been seen."""
        if db_seen_sequences is None:
            return False
        return seq_hash in db_seen_sequences

    def _write_fasta_record(self, handle, name: str, seq: str) -> None:
        """Write a FASTA record to a gzip file handle."""
        handle.write(f">{name}\n{seq}\n".encode())

    def _get_header(self, record) -> str:
        """Get the full header including comment from a sequence record."""
        if record.comment is None:
            header = record.name
        else:
            header = "_".join([record.name, record.comment])

        return header.replace(" ", "_")

    def _find_matching_databases(
        self, marker: str, databases: dict[str, dict[str, str]]
    ) -> list[str]:
        """Find all databases whose marker matches the given marker."""
        return [
            db_name
            for db_name, marker_info in databases.items()
            if self._marker_matches(marker, marker_info["marker"])
        ]

    def process(self, fasta: Path, databases: dict[str, dict[str, str]]) -> dict:
        """Process a FASTA file, splitting by markers with optional deduplication."""
        markers = [info["marker"] for info in databases.values()]
        logger.info(f"Splitting {fasta.name} by markers: {', '.join(markers)}")

        # Divide threads among writers to avoid oversubscription
        n_databases = len(databases)
        threads_per_writer = (
            self._threads // n_databases if n_databases <= self._threads else 1
        )

        # Open output files for each database with parallel gzip compression
        output_handles = {
            db: pgzip.open(
                self._tmpdir / f"{db}.fa.gz",
                "wb",
                thread=threads_per_writer,
                blocksize=2 * 10**8,
            )
            for db in databases.keys()
        }

        # Track which sequences are seen if we are deduplicating
        seen_sequences: dict[str, set[str]] | None = (
            {db: set() for db in databases.keys()} if self._deduplicate else None
        )

        # Count records processed per database
        record_counts = {db: 0 for db in databases.keys()}

        try:
            seq_count = 0
            with pysam.FastxFile(str(fasta), persist=False) as fh:
                for record in fh:
                    seq_count += 1
                    if seq_count % 500000 == 0:
                        logger.info(f"Processed {seq_count:,} sequences")

                    # Skip records with no sequence
                    sequence = record.sequence
                    if sequence is None:
                        continue

                    # Early length check before any other processing
                    if len(sequence) < self._min_length:
                        continue

                    header = self._get_header(record)
                    marker, taxon, output_header = self.header_processor(header)

                    # Find matching databases first
                    matching_dbs = self._find_matching_databases(marker, databases)
                    if not matching_dbs:
                        continue

                    # Lazy hash computation - only compute if we have matches and are deduplicating
                    seq_hash: str | None = None
                    if self._deduplicate:
                        seq_hash = self._compute_sequence_hash(sequence)

                    # Write to appropriate database files
                    for db_name in matching_dbs:
                        db_seen = seen_sequences[db_name] if seen_sequences else None

                        # Check for duplicates
                        if seq_hash is not None and self._is_duplicate(
                            seq_hash, db_seen
                        ):
                            continue

                        # Track this sequence as seen
                        if (
                            self._deduplicate
                            and seen_sequences is not None
                            and seq_hash is not None
                        ):
                            seen_sequences[db_name].add(seq_hash)

                        record_counts[db_name] += 1
                        self._write_fasta_record(
                            output_handles[db_name], output_header, sequence
                        )

        except IOError as e:
            raise IOError(f"Error processing {fasta.name}: {e}") from e

        finally:
            for handle in output_handles.values():
                handle.close()

        # Log databases with no records
        for db, count in record_counts.items():
            if count == 0:
                logger.warning(f"No records found for {db}")

        # Build output dictionary with only databases that have records
        return {
            db: {
                **databases[db],
                "processed_fasta": (self._tmpdir / f"{db}.fa.gz").resolve(),
                "n_seqs": record_counts[db],
            }
            for db in databases.keys()
            if record_counts[db] > 0
        }
