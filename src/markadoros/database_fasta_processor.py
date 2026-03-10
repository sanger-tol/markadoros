import gzip
import hashlib
from pathlib import Path
from typing import Callable

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
    ):
        self._deduplicate = deduplicate
        self._min_length = min_length
        self._tmpdir = tmpdir or Path.cwd()

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

    def _compute_sequence_hash(self, sequence: str | None) -> str:
        """Compute MD5 hash of canonical sequence."""
        if sequence is None:
            raise ValueError("Error: Input sequence is null!")

        canonical_seq = self._get_canonical_sequence(sequence)
        return hashlib.md5(canonical_seq.encode()).hexdigest()

    def _should_include_record(
        self,
        marker: str,
        sequence: str | None,
        seq_hash: str | None,
        db_marker: str,
        db_seen_sequences: set[str] | None,
    ) -> bool:
        """Check if a record should be included in the database."""
        if sequence is None:
            raise ValueError("Error: Input sequence is null!")

        # Check if marker matches
        if not marker.startswith(db_marker) and marker != db_marker:
            return False

        # Check minimum length
        if len(sequence) < self._min_length:
            return False

        # Check for duplicates
        if seq_hash is not None and db_seen_sequences is not None:
            if self._deduplicate and seq_hash in db_seen_sequences:
                return False

        return True

    def _write_fasta_record(self, handle, name: str, seq: str | None) -> None:
        """Write a FASTA record to a gzip file handle."""
        if seq is None:
            raise ValueError("Error: Input sequence is null!")

        handle.write(b">")
        handle.write(name.encode() + b"\n")
        handle.write(seq.encode() + b"\n")

    def _get_header(self, record) -> str:
        """Get the full header including comment from a sequence record."""
        if record.comment is None:
            header = record.name
        else:
            header = "_".join([record.name, record.comment])

        return header.replace(" ", "_")

    def process(self, fasta: Path, databases: dict[str, dict[str, str]]) -> dict:
        """Process a FASTA file, splitting by markers with optional deduplication."""
        markers = [info["marker"] for info in databases.values()]
        logger.info(f"Splitting {fasta.name} by markers: {', '.join(markers)}")

        # Open output files for each database
        output_handles = {
            db: gzip.open(self._tmpdir / f"{db}.fa.gz", "wb") for db in databases.keys()
        }

        # Track which sequences are seen if we are deduplicating
        seen_sequences = (
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

                    header = self._get_header(record)
                    marker, taxon, output_header = self.header_processor(header)

                    # Compute sequence hash for deduplication
                    seq_hash = None
                    if self._deduplicate:
                        seq_hash = self._compute_sequence_hash(record.sequence)

                    # Write to appropriate database files
                    for db_name, marker_info in databases.items():
                        if not self._should_include_record(
                            marker,
                            record.sequence,
                            seq_hash,
                            marker_info["marker"],
                            seen_sequences[db_name]
                            if seen_sequences is not None
                            else None,
                        ):
                            continue

                        if self._deduplicate and seen_sequences is not None:
                            seen_sequences[db_name].add(seq_hash)

                        record_counts[db_name] += 1
                        self._write_fasta_record(
                            output_handles[db_name], output_header, record.sequence
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
