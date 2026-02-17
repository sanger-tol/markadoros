import gzip
import io
import shutil
from contextlib import redirect_stdout
from pathlib import Path

import click
import pysam
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig
from pymmseqs.config import SearchConfig as MMSeqsSearchConfig
from pymmseqs.parsers import SearchParser as MMSeqsSearchParser

from markadoros.spades import SpadesRunner


class ShortReadAnalyser:
    def __init__(
        self,
        outdir: Path,
        prefix: str,
        threads: int,
        assembler: SpadesRunner,
        rmtmp: bool = False,
    ):
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True)

        ## Tempdir staging
        self.rmtmp = False
        self.tmpdir = self.outdir / "tmp"
        self.tmpdir.mkdir(parents=True, exist_ok=True)
        self.rmtmp = rmtmp

        ## Parameters
        self.threads = threads
        self.assembler = assembler
        self.prefix = prefix

    def _get_marker_from_db(self, db: Path) -> str:
        """Extract marker name from database path."""
        return db.parent.name

    def _extract_reads(
        self, search_result: MMSeqsSearchParser, input_reads: Path, output_path: Path
    ) -> Path:
        """Extract aligned reads to file. Returns path to output file."""
        with redirect_stdout(io.StringIO()):
            aligned_reads = {read["query"] for read in search_result.to_gen()}

        reads_extracted = 0
        with (
            pysam.FastxFile(input_reads) as fin,
            gzip.open(output_path, mode="wt") as fout,
        ):
            for read in fin:
                reads_extracted += 1
                if read.name in aligned_reads:
                    fout.write(str(read) + "\n")

        if reads_extracted == 0:
            return None

        return output_path

    def _create_db(self, fasta: Path, outdb: Path):
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
        alignment_db: str,
        tmp_dir: Path,
        min_seq_id: float = 0.9,
        min_aln_len: float = 50,
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
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
            )
            search.run()

        return MMSeqsSearchParser(search)

    def _assemble_reads(self, reads: Path, spades_outdir) -> Path | None:
        """Assemble reads and create database. Returns path to assembly DB."""

        contigs = self.assembler.assemble(reads, spades_outdir)

        if not contigs:
            return None

        return contigs

    def analyse_short_reads(
        self,
        input_reads: Path,
        marker: str,
        db: Path,
        min_seq_id: float = 0.96,
        min_aln_len: int = 100,
    ):
        """Run complete analysis pipeline."""
        # Create input MMSeqs DB
        input_db = self._create_db(input_reads, self.tmpdir / "input" / "db")

        # Search reads against DB
        click.echo(
            f"Searching reads with MMseqs2 against database: {marker}... ", nl=False
        )
        reads_search = self._search(
            query_db=input_db,
            target_db=db,
            alignment_db=self.tmpdir / marker / "search_reads" / "db",
            tmp_dir=(self.tmpdir / marker / "mmseqs_tmp" / "db"),
        )
        click.echo("done")

        # Extract aligned reads
        click.echo("Extracting aligned reads... ", nl=False)
        aligned_reads = self._extract_reads(
            reads_search,
            input_reads,
            self.tmpdir / marker / f"{Path(input_reads).stem}.match.fq.gz",
        )

        if aligned_reads is None:
            click.echo(f"No reads aligned to database {db}", err=True)
            return None

        click.echo("done")

        # Assemble reads and search contigs
        click.echo("Assembling extracted reads with SPAdes... ", nl=False)
        contigs = self._assemble_reads(
            aligned_reads,
            self.tmpdir / marker / "spades.out",
        )

        if not contigs or contigs.stat().st_size == 0:
            click.echo(f"No contigs assembled for database {db}", err=True)
            return None

        shutil.copy2(contigs, self.outdir / f"{self.prefix}.{marker}.contigs.fasta")
        click.echo("done")

        assembly_db = self._create_db(contigs, self.tmpdir / marker / "spades" / "db")

        ## Search for markers in assembled contigs
        click.echo(
            f"Searching assembled SPAdes contigs against database: {marker}... ",
            nl=False,
        )
        with redirect_stdout(io.StringIO()):
            contigs_search = self._search(
                query_db=db,
                target_db=assembly_db,
                alignment_db=self.tmpdir / "search_contigs" / "db",
                tmp_dir=(self.tmpdir / marker / "mmseqs_tmp" / "db"),
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
            )
            result = contigs_search.to_pandas()

        result["coverage"] = (
            result["target"]
            .str.extract(r"cov_(\d+\.?\d*)", expand=False)
            .astype(float)
            .round(2)
        )

        result.sort_values(
            by=["coverage", "query", "fident", "alnlen"], ascending=False
        ).to_csv(
            self.outdir / f"{self.prefix}.{marker}.result.tsv", index=False, sep="\t"
        )
        click.echo("done")

        return result
