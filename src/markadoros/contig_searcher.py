import io
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd
from loguru import logger
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig
from pymmseqs.config import SearchConfig as MMSeqsSearchConfig
from pymmseqs.parsers import SearchParser as MMSeqsSearchParser


class ContigSearcher:
    """Handles searching assembled contigs against marker databases."""

    def __init__(
        self,
        tmpdir: Path,
        threads: int,
    ):
        self.tmpdir = Path(tmpdir)
        self.threads = threads

    def search_contigs(
        self,
        contigs: Path,
        marker: str,
        marker_db: Path,
        min_seq_id: float = 0.96,
        min_aln_len: int = 400,
    ) -> pd.DataFrame | None:
        """Search contigs against a marker database.

        Args:
            contigs: Path to FASTA file with contigs
            marker: Marker/gene name for logging and temp dirs
            marker_db: Path to target marker database
            min_seq_id: Minimum sequence identity threshold
            min_aln_len: Minimum alignment length threshold

        Returns:
            DataFrame with search results, or None if search failed
        """
        logger.info(f"Searching contigs against {marker}...")

        # Create database from contigs
        contigs_db = self.tmpdir / marker / "contigs" / "db"
        contigs_db.parent.mkdir(parents=True, exist_ok=True)

        with redirect_stdout(io.StringIO()):
            db_config = CreateMMSeqsDBConfig(
                fasta_file=contigs.resolve(),
                sequence_db=contigs_db.resolve(),
                dbtype=2,
                v=3,
            )
            db_config.run()

        # Search contigs against target database
        alignment_db = self.tmpdir / marker / "search_contigs" / "db"
        alignment_db.parent.mkdir(parents=True, exist_ok=True)
        tmp_dir = self.tmpdir / marker / "mmseqs_tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        with redirect_stdout(io.StringIO()):
            search = MMSeqsSearchConfig(
                query_db=marker_db.resolve(),
                target_db=contigs_db.resolve(),
                alignment_db=alignment_db.resolve(),
                tmp_dir=tmp_dir.resolve(),
                search_type=3,
                v=3,
                threads=self.threads,
                max_accept=10,
                alignment_mode=2,
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
                strand=2,
                cov_mode=2,
                exact_kmer_matching=True,
            )
            search.run()

        search_result = MMSeqsSearchParser(search)
        with redirect_stdout(io.StringIO()):
            result_df = search_result.to_pandas()

        return result_df if not result_df.empty else None
