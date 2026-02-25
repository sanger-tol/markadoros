import gzip
import io
import shutil
from contextlib import redirect_stdout
from pathlib import Path

import click
import pandas as pd
import pysam
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig
from pymmseqs.config import SearchConfig as MMSeqsSearchConfig
from pymmseqs.parsers import SearchParser as MMSeqsSearchParser

from markadoros.assembler_runners import AssemblerRunner


class ReadAnalyser:
    def __init__(
        self,
        outdir: Path,
        prefix: str,
        threads: int,
        assembler: AssemblerRunner,
    ):
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True)

        ## Tempdir staging
        self.tmpdir = self.outdir / "tmp"
        self.tmpdir.mkdir(parents=True, exist_ok=True)

        ## Parameters
        self.threads = threads
        self.assembler = assembler
        self.prefix = prefix

    def _extract_reads(
        self, search_result: MMSeqsSearchParser, input_reads: Path, output_path: Path
    ) -> Path | None:
        """Extract aligned reads to file. Returns path to output file."""
        with redirect_stdout(io.StringIO()):
            aligned_reads = {read["query"] for read in search_result.to_gen()}

        reads_extracted = 0
        with (
            pysam.FastxFile(str(input_reads)) as fin,
            gzip.open(output_path, mode="wt") as fout,
        ):
            for read in fin:
                if read.name in aligned_reads:
                    reads_extracted += 1
                    fout.write(str(read) + "\n")

        if reads_extracted == 0:
            return None

        return output_path

    def _create_db(self, fasta: Path, outdb: Path) -> Path:
        """
        Create an MMSeqs2 DB.
        """
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
        s: float = 5.7,
        max_seqs: int = 300,
        alignment_mode: int = 2,
        min_seq_id: float = 0,
        min_aln_len: int = 0,
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
                s=s,
                max_seqs=max_seqs,
                alignment_mode=alignment_mode,
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
            )
            search.run()

        return MMSeqsSearchParser(search)

    def _filter_reads(self, input_reads: Path, marker: str, db: Path) -> Path | None:
        """Pre-filter reads against database."""
        input_db = self._create_db(input_reads, self.tmpdir / "input" / "db")

        # Run a loose pre-filtering search of reads against DB
        click.echo(f"Searching reads with MMseqs2 against database: {marker}... ")
        reads_search = self._search(
            query_db=input_db,
            target_db=db,
            alignment_db=self.tmpdir / marker / "search_reads" / "db",
            tmp_dir=self.tmpdir / marker / "mmseqs_tmp",
            s=1,
            max_seqs=10,
            alignment_mode=1,
        )

        # Extract aligned reads
        click.echo("Extracting aligned reads... ")
        aligned_reads = self._extract_reads(
            reads_search,
            input_reads,
            self.tmpdir / marker / f"{Path(input_reads).stem}.match.fq.gz",
        )

        return aligned_reads

    def _assemble_reads(self, aligned_reads: Path, marker: str) -> Path | None:
        """Assemble filtered reads into contigs"""
        click.echo("Assembling extracted reads with SPAdes... ")
        contigs = self.assembler.assemble(
            aligned_reads,
            self.tmpdir / marker / "spades.out",
        )

        if not contigs or contigs.stat().st_size == 0:
            return None

        shutil.copy2(contigs, self.outdir / f"{self.prefix}.{marker}.contigs.fasta")
        assembly_db = self._create_db(contigs, self.tmpdir / marker / "spades" / "db")

        return assembly_db

    def _search_contigs(
        self,
        contigs_db: Path,
        marker: str,
        db: Path,
        min_seq_id: float = 0.96,
        min_aln_len: int = 400,
    ) -> pd.DataFrame:
        click.echo(f"Searching assembled SPAdes contigs against database: {marker}... ")
        with redirect_stdout(io.StringIO()):
            contigs_search = self._search(
                query_db=db,
                target_db=contigs_db,
                alignment_db=self.tmpdir / marker / "search_contigs" / "db",
                tmp_dir=(self.tmpdir / marker / "mmseqs_tmp" / "db"),
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
            )
            result = contigs_search.to_pandas()

        return result

    def _save_result(self, result: pd.DataFrame, marker: str) -> pd.DataFrame:
        result["coverage"] = (
            result["target"]
            .str.extract(r"cov_(\d+\.?\d*)", expand=False)
            .astype(float)
            .round(2)
        )

        result.sort_values(by=["coverage", "fident", "alnlen"], ascending=False).to_csv(
            self.outdir / f"{self.prefix}.{marker}.result.tsv", index=False, sep="\t"
        )

        return result

    def analyse_short_reads(
        self,
        input_reads: Path,
        marker: str,
        db: Path,
        min_seq_id: float = 0.96,
        min_aln_len: int = 400,
    ) -> pd.DataFrame | None:
        """Run complete analysis pipeline."""
        aligned_reads = self._filter_reads(input_reads, marker, db)

        if aligned_reads is None:
            click.echo(f"No reads aligned to database {db}", err=True)
            return None

        assembled_contigs = self._assemble_reads(aligned_reads, marker)

        if assembled_contigs is None:
            click.echo(f"No contigs assembled for database {db}", err=True)
            return None

        contig_search_df = self._search_contigs(
            assembled_contigs,
            marker,
            db,
            min_seq_id=min_seq_id,
            min_aln_len=min_aln_len,
        )
        result = self._save_result(contig_search_df, marker)

        return result
