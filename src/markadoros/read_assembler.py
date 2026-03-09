import io
from contextlib import redirect_stdout
from pathlib import Path

import bgzip
import pysam
from loguru import logger
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

    def _extract_reads(
        self, search_result: MMSeqsSearchParser, input_reads: Path, output_path: Path
    ) -> tuple[int, Path | None]:
        """Extract aligned reads to file. Returns path to output file or None if no reads."""
        with redirect_stdout(io.StringIO()):
            aligned_reads = {read["query"] for read in search_result.to_gen()}

        reads_extracted = 0
        with (
            pysam.FastxFile(str(input_reads)) as fin,
            open(output_path, mode="wb") as fout,
        ):
            with bgzip.BGZipWriter(fout, num_threads=self.threads) as writer:
                for read in fin:
                    if read.name in aligned_reads:
                        reads_extracted += 1
                        writer.write((str(read) + "\n").encode("utf-8"))

        if reads_extracted == 0:
            return 0, None

        return reads_extracted, output_path

    def _filter_reads(
        self, input_reads: Path, marker: str, db: Path
    ) -> tuple[int, Path | None]:
        """Pre-filter reads against database and extract aligned reads."""
        input_db_path = self.tmpdir / "input" / "db"
        input_db_path.parent.mkdir(parents=True, exist_ok=True)

        with redirect_stdout(io.StringIO()):
            db_config = CreateMMSeqsDBConfig(
                fasta_file=input_reads.resolve(),
                sequence_db=input_db_path.resolve(),
                dbtype=2,
                v=3,
            )
            db_config.run()

        search_result_db = self.tmpdir / marker / "search_reads" / "db"
        search_result_db.parent.mkdir(parents=True, exist_ok=True)
        tmp_dir = self.tmpdir / marker / "mmseqs_tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        with redirect_stdout(io.StringIO()):
            search_config = MMSeqsSearchConfig(
                query_db=input_db_path.resolve(),
                target_db=db.resolve(),
                alignment_db=search_result_db.resolve(),
                tmp_dir=tmp_dir.resolve(),
                search_type=3,
                v=3,
                threads=self.threads,
                s=1,
                max_seqs=10,
                max_accept=1,
                alignment_mode=1,
                min_seq_id=0,
                min_aln_len=0,
                exact_kmer_matching=True,
            )
            search_config.run()

        reads_search = MMSeqsSearchParser(search_config)

        n_extracted_reads, aligned_reads = self._extract_reads(
            reads_search,
            input_reads,
            self.tmpdir / marker / f"{Path(input_reads).stem}.match.fq.gz",
        )

        return n_extracted_reads, aligned_reads

    def _assemble_reads(self, aligned_reads: Path, marker: str) -> Path | None:
        """Assemble filtered reads into contigs."""
        contigs = self.assembler.assemble(
            aligned_reads,
            self.tmpdir / marker / "assembly.out",
        )

        if not contigs or contigs.stat().st_size == 0:
            return None

        return contigs

    def assemble(
        self,
        input_reads: Path,
        marker: str,
        db: Path,
    ) -> tuple[int, Path | None]:
        """Filter reads and assemble them into contigs.

        Args:
            input_reads: Path to input reads file
            marker: Marker/gene name
            db: Path to target database for pre-filtering

        Returns:
            Path to assembled contigs FASTA file, or None if assembly failed
        """
        # Filter reads against target database
        n_aligned_reads, aligned_reads = self._filter_reads(input_reads, marker, db)

        logger.info(f"Searching reads against {marker}...")
        if aligned_reads is None:
            logger.error(f"No reads aligned to {marker}")
            return n_aligned_reads, None
        else:
            logger.info(f"Extracted {n_aligned_reads} reads aligning to {marker}!")

        # Assemble filtered reads
        logger.info(f"Assembling reads for {marker}...")
        assembled_contigs = self._assemble_reads(aligned_reads, marker)

        if assembled_contigs is None:
            logger.error(f"Failed to assemble contigs for {marker}")
            return n_aligned_reads, None

        return n_aligned_reads, assembled_contigs
