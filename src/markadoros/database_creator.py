import shutil
from pathlib import Path
from typing import Callable

from loguru import logger

from markadoros.database_fasta_processor import DatabaseFASTAProcessor
from markadoros.database_index import DatabaseIndex
from markadoros.mmseqs_database_builder import MMSeqsDatabaseBuilder


class DatabaseCreator:
    """Orchestrates the creation of marker databases from FASTA files."""

    def __init__(
        self,
        outdir: Path,
        header_processor: Callable[[str], tuple[str, str, str]] | None = None,
        deduplicate: bool = True,
        cluster: bool = False,
        cluster_min_seq_id: float = 0.99,
        cluster_coverage: float = 0.8,
        min_length: int = 200,
        create_index: bool = False,
        threads: int = 1,
    ) -> None:
        # Stage output directory
        self._outdir = Path(outdir)
        if not self._outdir.exists():
            self._outdir.mkdir(parents=True, exist_ok=True)

        # Setup temporary directory
        self._tmpdir = self._outdir / "tmp"
        self._tmpdir.mkdir(parents=True, exist_ok=True)

        self.deduplicate = deduplicate
        self.cluster = cluster
        self.min_length = min_length
        self.threads = threads

        # Initialize FastaProcessor
        self._fasta_processor = DatabaseFASTAProcessor(
            deduplicate=deduplicate,
            header_processor=header_processor,
            tmpdir=self._tmpdir,
            min_length=min_length,
            threads=threads,
        )

        # Initialize MarkerDatabaseBuilder and DatabaseIndex
        self._db_builder = MMSeqsDatabaseBuilder(
            self._outdir,
            self._tmpdir,
            threads=threads,
            cluster=cluster,
            cluster_perc_id=cluster_min_seq_id,
            cluster_coverage=cluster_coverage,
            create_index=create_index,
        )
        self._db_index = DatabaseIndex(self._outdir / "db.json")

    def create_marker_database(
        self,
        fasta: Path,
        markers: list[str],
        prefix: str,
    ) -> None:
        """Create MMSeqs2 databases from a FASTA file."""
        # Process FASTA file
        processed_dict = self._fasta_processor.process(
            fasta=fasta,
            databases={f"{prefix}.{x}": {"marker": x} for x in markers},
        )

        # Build MMSeqs database for each non-empty FASTA output processed
        for database, params in processed_dict.items():
            params = self._db_builder.build(database, params)

            # Add the database path and build FASTA to the db dict, and
            # remove the temporary FASTA file from it
            db_entry = {
                **params,
                "built_from": str(fasta.resolve()),
                "deduplicated": self.deduplicate,
                "clustered": self.cluster,
            }

            # Add the database entry to the index
            self._db_index.add_entry(database, db_entry)

        logger.info(f"Created {len(processed_dict)} databases! Exiting.")

    def cleanup(self) -> None:
        """Remove temporary files."""
        shutil.rmtree(self._tmpdir)
