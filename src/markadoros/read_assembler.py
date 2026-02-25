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

from markadoros.assembler_runners import AssemblerRunner


class ReadAssembler:
    """Handles read filtering and assembly into contigs."""

    def __init__(
        self,
        outdir: Path,
        tmpdir: Path,
        threads: int,
        prefix: str,
        assembler: AssemblerRunner,
    ):
        self.outdir = Path(outdir)
        self.tmpdir = Path(tmpdir)
        self.threads = threads
        self.prefix = prefix
        self.assembler = assembler

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
        s: float = 5.7,
        max_seqs: int = 300,
        max_accept: int = 100,
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
                max_accept=max_accept,
                alignment_mode=alignment_mode,
                min_seq_id=min_seq_id,
                min_aln_len=min_aln_len,
            )
            search.run()

        return MMSeqsSearchParser(search)

    def _extract_reads(
        self, search_result: MMSeqsSearchParser, input_reads: Path, output_path: Path
    ) -> Path | None:
        """Extract aligned reads to file. Returns path to output file or None if no reads."""
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

    def _filter_reads(self, input_reads: Path, marker: str, db: Path) -> Path | None:
        """Pre-filter reads against database and extract aligned reads."""
        input_db = self._create_db(input_reads, self.tmpdir / "input" / "db")

        click.echo(f"Searching reads against {marker}...")
        reads_search = self._search(
            query_db=input_db,
            target_db=db,
            alignment_db=self.tmpdir / marker / "search_reads" / "db",
            tmp_dir=self.tmpdir / marker / "mmseqs_tmp",
            s=1,
            max_seqs=10,
            max_accept=1,
            alignment_mode=1,
        )

        click.echo("Extracting aligned reads...")
        aligned_reads = self._extract_reads(
            reads_search,
            input_reads,
            self.tmpdir / marker / f"{Path(input_reads).stem}.match.fq.gz",
        )

        return aligned_reads

    def _assemble_reads(self, aligned_reads: Path, marker: str) -> Path | None:
        """Assemble filtered reads into contigs."""
        click.echo(f"Assembling reads for {marker}...")
        contigs = self.assembler.assemble(
            aligned_reads,
            self.tmpdir / marker / "assembly.out",
        )

        if not contigs or contigs.stat().st_size == 0:
            return None

        # Copy contigs to output directory
        output_contigs = self.outdir / f"{self.prefix}.{marker}.contigs.fasta"
        shutil.copy2(contigs, output_contigs)

        return output_contigs

    def assemble(
        self,
        input_reads: Path,
        marker: str,
        db: Path,
    ) -> Path | None:
        """Filter reads and assemble them into contigs.

        Args:
            input_reads: Path to input reads file
            marker: Marker/gene name
            db: Path to target database for pre-filtering

        Returns:
            Path to assembled contigs FASTA file, or None if assembly failed
        """
        # Filter reads against target database
        aligned_reads = self._filter_reads(input_reads, marker, db)

        if aligned_reads is None:
            click.echo(f"No reads aligned to {marker}", err=True)
            return None

        # Assemble filtered reads
        assembled_contigs = self._assemble_reads(aligned_reads, marker)

        if assembled_contigs is None:
            click.echo(f"Failed to assemble contigs for {marker}", err=True)
            return None

        return assembled_contigs
