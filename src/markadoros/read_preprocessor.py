import gzip
from pathlib import Path

import click
import pysam

from markadoros.utils import get_simple_name


class ReadPreprocessor:
    def __init__(
        self,
        input: Path,
        n_reads: int,
        outdir: Path,
    ):
        ## Input configuration
        self.input = input
        self.simple_name = get_simple_name(input)
        self.n_reads = n_reads
        self.subsampled_reads = None
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True)

    def _subsample_reads_cram(self) -> Path:
        """
        Read a CRAM file, and write out nreads reads to an interleaved FASTQ file.
        """
        outfile = self.outdir / f"{self.simple_name}.subsampled.fastq.gz"

        if not Path(str(self.input) + ".crai").exists():
            raise ValueError(f"CRAM file {self.input} is missing .crai index")

        cram = pysam.AlignmentFile(self.input, "rc", check_sq=False)
        with gzip.open(outfile, "wt") as f:
            for read in cram.head(self.n_reads):
                name_suffix = "/1" if read.is_read1 else "/2"
                f.write(f"@{read.query_name}{name_suffix}\n")
                f.write(f"{read.query_sequence}\n")
                f.write("+\n")
                f.write(f"{''.join(chr(q + 33) for q in read.query_qualities)}\n")

        cram.close()

        return outfile

    def _subsample_reads_fastx(self):
        """
        Read a FastQ file and write out the first N reads
        """
        outfile = self.outdir / f"{self.simple_name}.subsampled.fastq.gz"

        with pysam.FastxFile(self.input) as fin, gzip.open(outfile, mode="wt") as fout:
            for i, entry in enumerate(fin):
                if i >= self.n_reads:
                    break

                fout.write(str(entry) + "\n")

        return outfile

    def subsample_reads(self):
        """
        Get a subsample of reads and build an MMSeqs2 database from them
        """
        handlers = {
            ".cram": self._subsample_reads_cram,
            ".fq.gz": self._subsample_reads_fastx,
            ".fastq.gz": self._subsample_reads_fastx,
        }

        file_key = (
            "".join(self.input.suffixes[-2:])
            if self.input.suffix == ".gz"
            else self.input.suffix
        )

        if file_key not in handlers:
            raise ValueError(f"Unsupported input format: {self.input.name}")

        click.echo(f"Extracting the first {self.n_reads} reads from {self.input.name}")
        subsampled_reads = handlers[file_key]()

        return subsampled_reads
