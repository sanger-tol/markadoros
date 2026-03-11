import gzip
import io
import json
from contextlib import redirect_stdout
from pathlib import Path

import pysam
from loguru import logger
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig
from pymmseqs.config import CreateIndexConfig as MMSeqsCreateIndexConfig
from pymmseqs.config import EasyLinClustConfig as MMSeqsEasyLinClustConfig


class MMSeqsDatabaseBuilder:
    """Handles MMSeqs2 database creation and indexing."""

    def __init__(
        self,
        outdir: Path,
        tmpdir: Path,
        cluster: bool = False,
        threads: int = 1,
        cluster_perc_id: float = 0.99,
        cluster_coverage: float = 0.8,
    ):
        self._outdir = outdir
        self._tmpdir = tmpdir
        self.threads = threads
        self._cluster = cluster
        self._cluster_perc_id = cluster_perc_id
        self._cluster_coverage = cluster_coverage

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

    def _cluster_db(self, fasta_file: Path, cluster_prefix: Path) -> Path:
        """Cluster Fasta file."""
        with redirect_stdout(io.StringIO()):
            cluster_config = MMSeqsEasyLinClustConfig(
                fasta_files=fasta_file.resolve(),
                cluster_prefix=cluster_prefix.resolve(),
                tmp_dir=self._tmpdir.resolve(),
                threads=self.threads,
                min_seq_id=self._cluster_perc_id,
                c=self._cluster_coverage,
            )
            cluster_config.run()

        rep_seq_path = Path(str(cluster_prefix) + "_rep_seq.fasta").resolve()
        if not rep_seq_path.exists():
            raise FileNotFoundError(
                f"Expected clustered output not found: {rep_seq_path}"
            )

        return rep_seq_path

    def _index_db(self, db_path: Path) -> None:
        """Index MMSeqs2 database."""
        with redirect_stdout(io.StringIO()):
            index_config = MMSeqsCreateIndexConfig(
                sequence_db=db_path.resolve(),
                tmp_dir=self._tmpdir.resolve(),
                search_type=3,
                threads=self.threads,
            )
            index_config.run()

    def _count_taxa(self, db_path: Path) -> dict[str, int]:
        taxon_counts = {}
        with pysam.FastxFile(str(db_path)) as fasta:
            for record in fasta:
                if record.name is not None:
                    taxon = record.name.split("|")[2].replace("_", " ")
                else:
                    continue
                taxon_counts[taxon] = taxon_counts.get(taxon, 0) + 1
        return taxon_counts

    def _store_taxon_counts(
        self, taxon_index_path: Path, taxon_index: dict[str, int]
    ) -> None:
        with gzip.open(taxon_index_path, "wt") as f:
            json.dump(taxon_index, f, indent=2)

    def build(self, database: str, params: dict) -> tuple[Path, Path]:
        """Build and index an MMSeqs2 database."""
        db_path = self._outdir / database / "db"
        taxon_db_path = self._outdir / database / "taxon.json.gz"

        if self._cluster:
            logger.info(f"Clustering sequences for {database}...")
            cluster_prefix = self._tmpdir / (database + ".cluster")
            params["processed_fasta"] = self._cluster_db(
                fasta_file=params["processed_fasta"], cluster_prefix=cluster_prefix
            )

        logger.info(f"Building MMSeqs2 database for {database}... ")
        self._create_db(db_path, params["processed_fasta"])

        logger.info(f"Indexing MMSeqs2 database for {database}... ")
        self._index_db(db_path)

        logger.info(f"Counting taxa counts for {database}...")
        taxon_counts = self._count_taxa(params["processed_fasta"])
        self._store_taxon_counts(taxon_db_path, taxon_counts)

        return (db_path, taxon_db_path)
