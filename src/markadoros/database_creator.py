import shutil
from pathlib import Path
from typing import Callable

from loguru import logger

from markadoros.database_fasta_processor import DatabaseFASTAProcessor
from markadoros.database_index import DatabaseIndex
from markadoros.marker_database_builder import MarkerDatabaseBuilder


class DatabaseCreator:
    """Orchestrates the creation of marker databases from FASTA files."""

    def __init__(
        self,
        outdir: Path,
        db_dict: dict,
        header_processor: Callable[[str], tuple[str, str]] | None = None,
    ) -> None:
        # Stage output directory
        self._outdir = Path(outdir)
        if not self._outdir.exists():
            self._outdir.mkdir(parents=True, exist_ok=True)

        # Setup temporary directory
        self._tmpdir = self._outdir / "tmp"
        self._tmpdir.mkdir(parents=True, exist_ok=True)

        # Extract database parameters
        required_keys = {"parameters", "databases"}
        missing_keys = required_keys - set(db_dict.keys())
        if missing_keys:
            raise KeyError(
                f"Missing required keys in database configuration: {missing_keys}"
            )

        parameters = db_dict["parameters"]
        databases = db_dict["databases"]

        self.deduplicated = parameters.get("deduplicate", False)

        # Initialize FastaProcessor
        self._fasta_processor = DatabaseFASTAProcessor(
            databases=databases,
            deduplicate=parameters.get("deduplicate", False),
            header_processor=header_processor,
            tmpdir=self._tmpdir,
        )

        # Initialize MarkerDatabaseBuilder and DatabaseIndex
        self._db_builder = MarkerDatabaseBuilder(self._outdir, self._tmpdir)
        self._db_index = DatabaseIndex(self._outdir / "db.json")

    def create_marker_database(self, fasta: Path) -> None:
        """Create MMSeqs2 databases from a FASTA file."""
        # Process FASTA file
        processed_dict = self._fasta_processor.process(fasta)

        # Build MMSeqs database for each non-empty FASTA output processed
        for database, params in processed_dict.items():
            db_path = self._db_builder.build(database, params)

            # Add the database path and build FASTA to the db dict, and
            # remove the temporary FASTA file from it
            db_entry = {
                **params,
                "db": str(db_path.resolve()),
                "built_from": str(fasta.resolve()),
                "deduplicated": self.deduplicated,
            }
            db_entry.pop("processed_fasta", None)

            # Add the database entry to the index
            self._db_index.add_entry(database, db_entry)

        logger.info(f"Created {len(processed_dict)} databases! Exiting.")

    def cleanup(self) -> None:
        """Remove temporary files."""
        shutil.rmtree(self._tmpdir)
