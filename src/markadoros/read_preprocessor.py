from pathlib import Path

import pysam
from isal import igzip
from loguru import logger

from markadoros.utils import get_simple_name


class ReadPreprocessor:
    def __init__(
        self,
        outdir: Path,
        threads: int = 1,
    ):
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True)

        self.threads = threads

    def _preprocess_reads_sam(
        self, input: Path, nreads: int | None = None, mode: str = "cram"
    ) -> Path:
        """
        Read a BAM/CRAM file, and write out nreads reads to an interleaved FASTQ file.
        """
        outfile = self.outdir / f"{get_simple_name(input)}.subsampled.fastq.gz"
        read_mode = "rc" if mode == "cram" else "rb"
        sam = pysam.AlignmentFile(
            str(input), read_mode, check_sq=False, require_index=True
        )

        with igzip.open(outfile, "wb") as fout:
            collected = []
            for count, read in enumerate(sam.fetch(".")):
                if nreads is not None and count >= nreads:
                    break
                name_suffix = "/1" if read.is_read1 else "/2"
                quals = (
                    bytes(q + 33 for q in read.query_qualities)
                    if read.query_qualities is not None
                    else b"+"
                )
                collected.append(
                    f"@{read.query_name}{name_suffix}\n{read.query_sequence}\n+\n".encode(
                        "utf-8"
                    )
                    + quals
                    + b"\n"
                )
                if len(collected) >= 100_000:
                    fout.writelines(collected)
                    collected.clear()
            if collected:
                fout.writelines(collected)

        sam.close()
        return outfile

    def _preprocess_reads_fastx(self, input: Path, nreads: int | None = None):
        is_fastq = input.suffix in (".fastq", ".fq") or str(input).endswith(".fastq.gz")
        out_ext = "fastq" if is_fastq else "fasta"
        outfile = self.outdir / f"{get_simple_name(input)}.subsampled.{out_ext}.gz"

        in_opener = igzip.open if str(input).endswith(".gz") else open

        with in_opener(input, "rb") as fin, igzip.open(outfile, "wb") as fout:
            if is_fastq:
                collected = []
                lines_needed = nreads * 4 if nreads is not None else float("inf")
                for line in fin:
                    collected.append(line)
                    if len(collected) >= lines_needed:
                        break
                fout.writelines(collected)
            else:
                count = 0
                collected = []
                for line in fin:
                    if line[0:1] == b">":
                        if nreads is not None and count == nreads:
                            break
                        count += 1
                    collected.append(line)
                fout.writelines(collected)

        return outfile

    def preprocess_reads(
        self,
        input_file: Path,
        n_reads: int | None = None,
    ) -> Path:
        """
        Get a subsample of reads and build an MMSeqs2 database from them
        """
        if not input_file.exists():
            raise FileNotFoundError(f"Input file {input_file} does not exist!")

        # Determine file type
        file_key = (
            "".join(input_file.suffixes[-2:])
            if input_file.suffix == ".gz"
            else input_file.suffix
        )

        # If no subsampling requested, return input directly for non-CRAM files
        if n_reads is None and file_key != ".cram":
            return input_file

        # Print appropriate status message
        action = "Converting" if file_key == ".cram" else "Extracting"
        read_limit = "all" if n_reads is None else f"first {n_reads}"
        logger.info(f"{action} {read_limit} reads from {input_file.name}")

        if file_key == ".cram" or file_key == ".bam":
            return self._preprocess_reads_sam(input_file, n_reads)
        else:
            return self._preprocess_reads_fastx(input_file, n_reads)
