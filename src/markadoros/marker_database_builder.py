import io
from contextlib import redirect_stdout
from pathlib import Path

import click
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

    def build(self, database: str, params: dict, fasta: Path) -> Path:
        """Build and index an MMSeqs2 database."""
        db_path = self._outdir / database / "db"

        click.echo(f"Building MMSeqs2 database for {database}... ")
        self._create_db(db_path, params.get("processed_fasta"))

        click.echo(f"Indexing MMSeqs2 database for {database}... ")
        self._index_db(db_path)

        return db_path
