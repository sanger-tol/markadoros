import gzip
from pathlib import Path

import click
import pysam

from markadoros.utils import get_simple_name


class ReadPreprocessor:
    def __init__(
        self,
        outdir: Path,
    ):
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True)

    def _preprocess_reads_cram(self, input: Path, nreads: int | None = None) -> Path:
        """
        Read a CRAM file, and write out nreads reads to an interleaved FASTQ file.
        """
        count = 0
        outfile = self.outdir / f"{get_simple_name(input)}.subsampled.fastq.gz"

        cram = pysam.AlignmentFile(str(input), "rc", check_sq=False, require_index=True)
        with gzip.open(outfile, "wt") as f:
            for read in cram.fetch("."):
                count += 1
                if nreads is not None and count > nreads:
                    break

                name_suffix = "/1" if read.is_read1 else "/2"
                f.write(f"@{read.query_name}{name_suffix}\n")
                f.write(f"{read.query_sequence}\n")
                f.write("+\n")
                if read.query_qualities is not None:
                    f.write(f"{''.join(chr(q + 33) for q in read.query_qualities)}\n")
                else:
                    f.write("+\n")

        cram.close()

        return outfile

    def _preprocess_reads_fastx(self, input: Path, nreads: int | None = None):
        """
        Read a FastQ file and write out the first N reads
        """
        outfile = self.outdir / f"{get_simple_name(input)}.subsampled.fastq.gz"

        with pysam.FastxFile(str(input)) as fin, gzip.open(outfile, mode="wt") as fout:
            for i, entry in enumerate(fin):
                if nreads is not None and i >= nreads:
                    break

                fout.write(str(entry) + "\n")

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
        click.echo(f"{action} {read_limit} reads from {input_file.name}")

        if file_key == ".cram":
            return self._preprocess_reads_cram(input_file, n_reads)
        else:
            return self._preprocess_reads_fastx(input_file, n_reads)
