from pathlib import Path

import bgzip
import pysam
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

    def _preprocess_reads_cram(self, input: Path, nreads: int | None = None) -> Path:
        """
        Read a CRAM file, and write out nreads reads to an interleaved FASTQ file.
        """
        count = 0
        outfile = self.outdir / f"{get_simple_name(input)}.subsampled.fastq.gz"

        cram = pysam.AlignmentFile(str(input), "rc", check_sq=False, require_index=True)
        with open(outfile, "wb") as f:
            with bgzip.BGZipWriter(f, num_threads=self.threads) as writer:
                for read in cram.fetch("."):
                    count += 1
                    if nreads is not None and count > nreads:
                        break

                    name_suffix = "/1" if read.is_read1 else "/2"
                    writer.write(f"@{read.query_name}{name_suffix}\n".encode("utf-8"))
                    writer.write(f"{read.query_sequence}\n".encode("utf-8"))
                    writer.write("+\n".encode("utf-8"))
                    if read.query_qualities is not None:
                        writer.write(
                            f"{''.join(chr(q + 33) for q in read.query_qualities)}\n".encode(
                                "utf-8"
                            )
                        )
                    else:
                        writer.write("+\n".encode("utf-8"))

        cram.close()

        return outfile

    def _preprocess_reads_fastx(self, input: Path, nreads: int | None = None):
        """
        Read a FastQ file and write out the first N reads
        """
        outfile = self.outdir / f"{get_simple_name(input)}.subsampled.fastq.gz"

        with pysam.FastxFile(str(input)) as fin, open(outfile, mode="wb") as fout:
            with bgzip.BGZipWriter(fout, num_threads=self.threads) as writer:
                for i, entry in enumerate(fin):
                    if nreads is not None and i >= nreads:
                        break

                    writer.write((str(entry) + "\n").encode("utf-8"))

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

        if file_key == ".cram":
            return self._preprocess_reads_cram(input_file, n_reads)
        else:
            return self._preprocess_reads_fastx(input_file, n_reads)
