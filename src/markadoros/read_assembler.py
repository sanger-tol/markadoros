import io
from contextlib import redirect_stdout
from pathlib import Path

from isal import igzip
from loguru import logger
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig
from pymmseqs.config import SearchConfig as MMSeqsSearchConfig
from pymmseqs.parsers import SearchParser as MMSeqsSearchParser

from markadoros.assembler_runners import AssemblerRunner
from markadoros.input_types import InputType


class ReadAssembler:
    """Handles read filtering and assembly into contigs."""

    def __init__(
        self,
        outdir: Path,
        tmpdir: Path,
        threads: int,
        prefix: str,
        assembler: AssemblerRunner,
        input_type: InputType,
    ):
        self.outdir = Path(outdir)
        self.tmpdir = Path(tmpdir)
        self.threads = threads
        self.prefix = prefix
        self.assembler = assembler
        self.input_type = input_type

    def _extract_reads(
        self, search_result: MMSeqsSearchParser, input_reads: Path, output_path: Path
    ) -> tuple[int, Path | None]:
        """Extract aligned reads to file. Returns path to output file or None if no reads.

        Supports both FASTA (>) and FASTQ (@) formats. Format is detected by the first character.
        """
        with redirect_stdout(io.StringIO()):
            aligned_reads = set(search_result.to_pandas()["query"])

        reads_extracted = 0
        in_opener = igzip.open if str(input_reads).endswith(".gz") else open

        with in_opener(input_reads, "rb") as fin, igzip.open(output_path, "wb") as fout:
            # Detect format from first character
            first_char = fin.read(1)
            if not first_char:
                return 0, None

            is_fastq = first_char == b"@"

            # Seek back to beginning for processing
            fin.seek(0)
            collected = []

            if is_fastq:
                reads_extracted = self._extract_fastq_records(
                    fin, fout, aligned_reads, collected
                )
            else:
                reads_extracted = self._extract_fasta_records(
                    fin, fout, aligned_reads, collected
                )

        if reads_extracted == 0:
            return 0, None

        return reads_extracted, output_path

    def _extract_fasta_records(
        self, fin, fout, aligned_reads: set, collected: list
    ) -> int:
        """Extract FASTA records for reads in aligned_reads set."""
        reads_extracted = 0
        current_id = None
        write_current = False

        for line in fin:
            if line[0:1] == b">":
                current_id = line[1:].split()[0].decode()
                write_current = current_id in aligned_reads
                if write_current:
                    reads_extracted += 1

            if write_current:
                collected.append(line)
                if len(collected) >= 100_000:
                    fout.writelines(collected)
                    collected.clear()

        if collected:
            fout.writelines(collected)

        return reads_extracted

    def _extract_fastq_records(
        self, fin, fout, aligned_reads: set, collected: list
    ) -> int:
        """Extract FASTQ records for reads in aligned_reads set.

        FASTQ format: 4 lines per record
        - Line 1: header starting with @
        - Line 2: sequence
        - Line 3: separator (usually +)
        - Line 4: quality scores
        """
        reads_extracted = 0

        while True:
            # Read FASTQ record (4 lines)
            header = fin.readline()
            if not header:
                break

            seq = fin.readline()
            sep = fin.readline()
            qual = fin.readline()

            if not all([header, seq, sep, qual]):
                break

            # Extract read ID from header (remove @ and split on whitespace)
            read_id = header[1:].split()[0].decode()

            if read_id in aligned_reads:
                reads_extracted += 1
                collected.extend([header, seq, sep, qual])

                if len(collected) >= 100_000:
                    fout.writelines(collected)
                    collected.clear()

        if collected:
            fout.writelines(collected)

        return reads_extracted

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

        search_result_db = self.tmpdir / "search_reads" / "db"
        search_result_db.parent.mkdir(parents=True, exist_ok=True)
        tmp_dir = self.tmpdir / "mmseqs_tmp"
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
                max_seqs=50,
                max_accept=1,
                alignment_mode=1,
                min_seq_id=0.5,
                min_aln_len=450
                if (
                    self.input_type == InputType.PACBIO
                    or self.input_type == InputType.ONT
                )
                else 50,
                exact_kmer_matching=True,
            )
            search_config.run()

        reads_search = MMSeqsSearchParser(search_config)

        n_extracted_reads, aligned_reads = self._extract_reads(
            reads_search,
            input_reads,
            self.tmpdir / f"{Path(input_reads).stem}.match.fq.gz",
        )

        return n_extracted_reads, aligned_reads

    def _assemble_reads(self, aligned_reads: Path, marker: str) -> Path | None:
        """Assemble filtered reads into contigs."""
        contigs = self.assembler.assemble(
            aligned_reads,
            self.tmpdir / "assembly.out",
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
        logger.info(f"Searching reads against {marker}...")
        n_aligned_reads, aligned_reads = self._filter_reads(input_reads, marker, db)

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
