import re
import subprocess
from pathlib import Path

from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig


class SpadesRunner:
    def __init__(self, threads: int, rna: bool = False):
        self.rna = rna
        self.threads = threads

        ## Check if SPAdes is installed
        self._check_install()

    def _check_install(self, mmseqs_path: str) -> None:
        """Check if spades.py is available and extract version."""
        try:
            result = subprocess.run(
                ["spades.py", "--version"], capture_output=True, text=True, check=True
            )

            if "SPAdes" not in result.stdout:
                raise ValueError("MMSeqs2 output doesn't contain expected header")

            version = result.stdout.removeprefix("SPAdes genome assembler v").strip()

            if not re.match(r"^\d+\.\d+\.\d+$", version):
                raise ValueError("Could not find SPAdes version in output")

            self.spades_version = version

        except FileNotFoundError:
            raise ValueError("SPAdes is not installed!")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SPAdes failed: {e}")

    def _run_spades(self, reads_fq: Path, outdir: Path):
        """
        Run SPAdes to get an assembly from a set of input reads.
        """
        try:
            spades_process_call = [
                "spades.py",
                "--rna" if self.rna else None,
                "-t",
                self.threads,
                "-s",
                str({reads_fq.resolve()}),
                "-o",
                str(outdir),
            ]

            subprocess.run(
                [arg for arg in spades_process_call if arg],
                capture_output=True,
                text=True,
                check=True,
            )
        except FileNotFoundError:
            raise FileNotFoundError("SPAdes is not installed!")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SPAdes failed: {e}")

        spades_contigs = (
            "hard_filtered_transcripts.fasta" if self.rrna else "contigs.fasta"
        )

        return outdir / spades_contigs

    def assemble(self, reads_fq: Path, outdir: Path):
        """
        Assemble reads with SPAdes
        """
        assembly = self._run_spades(reads_fq, outdir)

        return assembly
