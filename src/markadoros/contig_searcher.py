import io
from contextlib import redirect_stdout
from pathlib import Path

import click
import pandas as pd
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

    def _create_db(self, fasta: Path, outdb: Path) -> Path:
        """Create an MMSeqs2 database from a FASTA file."""
        with redirect_stdout(io.StringIO()):
            db_config = CreateMMSeqsDBConfig(
                fasta_file=fasta.resolve(),
                sequence_db=outdb.resolve(),
                dbtype=2,
                v=1,
            )
            db_config.run()

        return outdb.resolve()

    def _search(
        self,
        query_db: Path,
        target_db: Path,
        alignment_db: Path,
        tmp_dir: Path,
        start_sens: float = 1,
        sens: float = 7,
        sens_steps: int = 3,
        max_seqs: int = 300,
        max_accept: int = 2 ^ 15,
        alignment_mode: int = 2,
        min_seq_id: float = 0.96,
        min_aln_len: int = 400,
    ) -> MMSeqsSearchParser:
        """Search sequences against database."""
        with redirect_stdout(io.StringIO()):
            search = MMSeqsSearchConfig(
                query_db=query_db.resolve(),
                target_db=target_db.resolve(),
                alignment_db=alignment_db.resolve(),
                tmp_dir=tmp_dir.resolve(),
                search_type=3,
                v=1,
                threads=self.threads,
                start_sens=start_sens,
                s=sens,
                sens_steps=sens_steps,
                max_seqs=max_seqs,
                max_accept=max_accept,
                alignment_mode=alignment_mode,
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
            )
            search.run()

        return MMSeqsSearchParser(search)

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
            target_db: Path to target marker database
            min_seq_id: Minimum sequence identity threshold
            min_aln_len: Minimum alignment length threshold

        Returns:
            DataFrame with search results, or None if search failed
        """
        click.echo(f"Searching contigs against {marker}...")

        try:
            # Create database from contigs
            contigs_db = self._create_db(
                contigs, self.tmpdir / marker / "contigs" / "db"
            )

            # Search contigs against target database
            search_result = self._search(
                query_db=marker_db,
                target_db=contigs_db,
                alignment_db=self.tmpdir / marker / "search_contigs" / "db",
                tmp_dir=self.tmpdir / marker / "mmseqs_tmp",
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
                max_accept=2,
                start_sens=1,
                sens=7,
                sens_steps=3,
            )

            with redirect_stdout(io.StringIO()):
                result_df = search_result.to_pandas()

            return result_df if not result_df.empty else None

        except Exception as e:
            click.echo(f"Error searching contigs: {e}", err=True)
            return None
