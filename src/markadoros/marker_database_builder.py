import gzip
import io
import json
from contextlib import redirect_stdout
from pathlib import Path

from loguru import logger
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig
from pymmseqs.config import CreateIndexConfig as MMSeqsCreateIndexConfig


class MarkerDatabaseBuilder:
    """Handles MMSeqs2 database creation and indexing."""

    def __init__(self, outdir: Path, tmpdir: Path):
        self._outdir = outdir
        self._tmpdir = tmpdir

    def _create_db(self, db_path: Path, fasta_file: Path) -> None:
        """Create MMSeqs2 database from FASTA file."""
        with redirect_stdout(io.StringIO()):
            db_config = CreateMMSeqsDBConfig(
                fasta_file=fasta_file,
                sequence_db=db_path.resolve(),
                dbtype=2,
                v=1,
            )
            db_config.run()

    def _index_db(self, db_path: Path) -> None:
        """Index MMSeqs2 database."""
        with redirect_stdout(io.StringIO()):
            index_config = MMSeqsCreateIndexConfig(
                sequence_db=db_path.resolve(),
                tmp_dir=self._tmpdir.resolve(),
                search_type=3,
            )
            index_config.run()

    def _store_taxon_counts(
        self, taxon_index_path: Path, taxon_index: dict[str, int]
    ) -> None:
        with gzip.open(taxon_index_path, "wt") as f:
            json.dump(taxon_index, f, indent=2)

    def build(self, database: str, params: dict) -> tuple[Path, Path]:
        """Build and index an MMSeqs2 database."""
        db_path = self._outdir / database / "db"
        taxon_db_path = self._outdir / database / "taxon.json.gz"

        logger.info(f"Building MMSeqs2 database for {database}... ")
        self._create_db(db_path, params["processed_fasta"])

        logger.info(f"Indexing MMSeqs2 database for {database}... ")
        self._index_db(db_path)

        logger.info(f"Writing taxon counts for {database}...")
        self._store_taxon_counts(taxon_db_path, params["taxon_counts"])

        return (db_path, taxon_db_path)
