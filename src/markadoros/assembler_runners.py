import re
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path

from loguru import logger


class AssemblerRunner(ABC):
    """Base class for genome/transcript assemblers."""

    def __init__(self, threads: int):
        self.threads = threads
        self._check_install()

    @abstractmethod
    def _check_install(self) -> None:
        """Check if assembler is installed and extract version."""
        pass

    @abstractmethod
    def _get_command(self, reads_fq: Path, outdir: Path) -> list:
        """Return the command to run the assembler."""
        pass

    @abstractmethod
    def _get_output_file(self, outdir: Path) -> Path:
        """Return the path to the output assembly file."""
        pass

    def _process_contigs(self, contigs: Path) -> Path:
        return contigs

    def _run_assembler(self, reads_fq: Path, outdir: Path):
        """
        Generic method to run an assembler with input reads.
        """
        if not outdir.exists():
            outdir.mkdir(parents=True)

        try:
            command = self._get_command(reads_fq, outdir)
            log_file = outdir / f"{self.__class__.__name__.lower()}.log"

            with open(log_file, "w") as log:
                subprocess.run(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=True,
                    cwd=outdir,
                )

        except FileNotFoundError:
            raise FileNotFoundError(f"{self.__class__.__name__} is not installed!")
        except subprocess.CalledProcessError as e:
            logger.error(
                f"{self.__class__.__name__} failed with exit code {e.returncode}. Check {log_file} for details."
            )
            return None

        output_file = self._get_output_file(outdir)
        if not output_file.exists():
            logger.error(
                f"{self.__class__.__name__} did not produce expected output. Check {log_file} for details."
            )
            return None

        processed_output = self._process_contigs(output_file)

        return processed_output

    def assemble(self, reads_fq: Path, outdir: Path):
        """
        Assemble reads with this assembler.
        """
        return self._run_assembler(reads_fq, outdir)


class SpadesRunner(AssemblerRunner):
    def __init__(self, threads: int, rna: bool = False):
        self.rna = rna
        super().__init__(threads)

    def _check_install(self) -> None:
        """Check if spades.py is available and extract version."""
        try:
            result = subprocess.run(
                ["spades.py", "--version"], capture_output=True, text=True, check=True
            )

            if "SPAdes" not in result.stdout:
                raise ValueError("SPAdes output doesn't contain expected header")
        except FileNotFoundError:
            raise ValueError("SPAdes is not installed!")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SPAdes failed: {e}")

    def _get_command(self, reads_fq: Path, outdir: Path) -> list:
        """Return the SPAdes command."""
        command = [
            "spades.py",
            "--rna" if self.rna else None,
            "-t",
            str(self.threads),
            "-s",
            str(reads_fq.resolve()),
            "-o",
            ".",
        ]
        return [arg for arg in command if arg]

    def _get_output_file(self, outdir: Path) -> Path:
        """Return the SPAdes output file path."""
        spades_contigs = (
            "hard_filtered_transcripts.fasta" if self.rna else "contigs.fasta"
        )
        return outdir / spades_contigs


class HifiasmRunner(AssemblerRunner):
    def __init__(self, threads: int, ont: bool = False):
        self.ont = ont
        super().__init__(threads)

    def _check_install(self) -> None:
        """Check if hifiasm is available and extract version."""
        try:
            result = subprocess.run(
                ["hifiasm", "--version"], capture_output=True, text=True, check=True
            )

            version = result.stdout.strip()
            self.hifiasm_version = version

        except FileNotFoundError:
            raise ValueError("hifiasm is not installed!")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"hifiasm failed: {e}")

    def _get_command(self, reads_fq: Path, outdir: Path) -> list:
        """Return the hifiasm command."""
        command = [
            "hifiasm",
            "-t",
            str(self.threads),
            "--ont" if self.ont else "",
            "-o",
            "asm",
            str(reads_fq.resolve()),
        ]

        return [arg for arg in command if arg]

    def _process_contigs(self, contigs: Path) -> Path:
        outfile = self._gfa_to_fasta(contigs)
        return outfile

    def _get_output_file(self, outdir: Path) -> Path:
        """Return the hifiasm output file path (primary assembly)."""
        return outdir / "asm.bp.p_ctg.gfa"

    def _gfa_to_fasta(self, gfa_file: Path) -> Path:
        """Convert GFA format to FASTA, extracting coverage from rd:i: field."""
        if not gfa_file.exists():
            raise FileNotFoundError(f"GFA file not found: {gfa_file}")

        outfile = gfa_file.parent / (gfa_file.stem + ".fa")

        with (
            open(gfa_file, "r") as fin,
            open(outfile, "w") as fout,
        ):
            for line in fin:
                if line.startswith("S"):
                    line_ = line.strip().split("\t")
                    header = f">{line_[1]}"
                    for field in line_:
                        match = re.search(r"rd:i:(\d+)", field)
                        if match:
                            header += f"_cov_{match.group(1)}"
                            break

                    fout.write(f"{header}\n{line_[2]}\n")

        return outfile
