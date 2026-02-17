import re
import subprocess
from pathlib import Path

import click


class SpadesRunner:
    def __init__(self, threads: int, rna: bool = False):
        self.rna = rna
        self.threads = threads

        ## Check if SPAdes is installed
        self._check_install()

    def _check_install(self) -> None:
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
        if not outdir.exists():
            outdir.mkdir(parents=True)

        try:
            spades_process_call = [
                "spades.py",
                "--rna" if self.rna else None,
                "-t",
                str(self.threads),
                "-s",
                str(reads_fq.resolve()),
                "-o",
                str(outdir),
            ]

            log_file = outdir / "spades.log"

            with open(log_file, "w") as log:
                subprocess.run(
                    [arg for arg in spades_process_call if arg],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=True,
                )

        except FileNotFoundError:
            raise FileNotFoundError("SPAdes is not installed!")
        except subprocess.CalledProcessError as e:
            click.echo(
                f"SPAdes failed with exit code {e.returncode}. Check {outdir}/spades.log for details.",
                err=True,
            )
            return None

        spades_contigs = (
            "hard_filtered_transcripts.fasta" if self.rna else "contigs.fasta"
        )

        output_file = outdir / spades_contigs
        if not output_file.exists():
            click.echo(
                f"SPAdes did not produce {spades_contigs}. Check {outdir}/spades.log for details.",
                err=True,
            )
            return None

        return outdir / spades_contigs

    def assemble(self, reads_fq: Path, outdir: Path):
        """
        Assemble reads with SPAdes
        """
        assembly = self._run_spades(reads_fq, outdir)

        return assembly
