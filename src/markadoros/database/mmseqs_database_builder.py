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
        create_index: bool = False,
    ):
        self._outdir = outdir
        self._tmpdir = tmpdir
        self.threads = threads
        self._cluster = cluster
        self._cluster_perc_id = cluster_perc_id
        self._cluster_coverage = cluster_coverage
        self._create_index = create_index

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

    def _cluster_db(
        self,
        fasta_file: Path,
        cluster_prefix: Path,
        cluster_perc_id: float,
        cluster_coverage: float,
    ) -> Path:
        """Cluster Fasta file."""
        with redirect_stdout(io.StringIO()):
            cluster_config = MMSeqsEasyLinClustConfig(
                fasta_files=fasta_file.resolve(),
                cluster_prefix=cluster_prefix.resolve(),
                tmp_dir=self._tmpdir.resolve(),
                threads=self.threads,
                min_seq_id=cluster_perc_id,
                c=cluster_coverage,
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

    def _count_records(self, db_path: Path) -> tuple[int, dict[str, int]]:
        record_counts = 0
        taxon_counts = {}
        with pysam.FastxFile(str(db_path)) as fasta:
            for record in fasta:
                if record.name is None:
                    continue
                elif record.comment is None:
                    header = record.name
                else:
                    header = " ".join([record.name, record.comment])

                record_counts += 1
                taxon = header.split("|")[2].replace("_", " ")
                taxon_counts[taxon] = taxon_counts.get(taxon, 0) + 1

        return record_counts, taxon_counts

    def _write_taxon_counts(
        self, taxon_index_path: Path, taxon_index: dict[str, int]
    ) -> None:
        with gzip.open(taxon_index_path, "wt") as f:
            json.dump(taxon_index, f, indent=2)

    def build(self, database: str, params: dict) -> dict[str, str]:
        """Build and index an MMSeqs2 database."""
        params["db"] = self._outdir / database / "db"
        params["taxon_db"] = self._outdir / database / "taxon.json.gz"

        if self._cluster:
            logger.info(f"Clustering sequences for {database}...")
            cluster_prefix = self._tmpdir / (database + ".cluster")
            db_input_fasta = self._cluster_db(
                fasta_file=params["processed_fasta"],
                cluster_prefix=cluster_prefix,
                cluster_perc_id=self._cluster_perc_id,
                cluster_coverage=self._cluster_coverage,
            )
        else:
            db_input_fasta = params["processed_fasta"]

        logger.info(f"Building MMSeqs2 database for {database}... ")
        self._create_db(params["db"], db_input_fasta)

        if self._create_index:
            logger.info(f"Indexing MMSeqs2 database for {database}... ")
            self._index_db(params["db"])

        logger.info(f"Counting records for {database}...")
        record_counts, taxon_counts = self._count_records(db_input_fasta)
        params["n_seqs"] = record_counts
        self._write_taxon_counts(params["taxon_db"], taxon_counts)

        params.pop("processed_fasta")

        return params
