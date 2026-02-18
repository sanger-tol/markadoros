import shutil
from pathlib import Path

import click

from markadoros.assembler import HifiasmRunner, SpadesRunner
from markadoros.read_analyser import ReadAnalyser
from markadoros.read_preprocessor import ReadPreprocessor


class SearchPipeline:
    """Orchestrates the barcode search workflow."""

    def __init__(
        self,
        outdir: Path,
        threads: int,
        database_index: dict,
        platform: str,
        cleanup: bool = True,
    ):
        self.outdir = Path(outdir)
        self.threads = threads
        self.database_index = database_index

        if platform == "illumina":
            self.assembler = SpadesRunner(threads=self.threads, rna=False)
        elif platform == "rnaseq" or platform == "illumina_rnaseq":
            self.assembler = SpadesRunner(threads=self.threads, rna=True)
        elif platform == "pacbio_hifi" or platform == "pb":
            self.assembler = HifiasmRunner(threads=self.threads, ont=False)
        elif platform == "oxford_nanopore" or platform == "ont":
            self.assembler = HifiasmRunner(threads=self.threads, ont=True)
        else:
            raise ValueError(f"Unsupported platform: {platform}")

    def _print_result_summary(self, result) -> None:
        """Print a summary of the top search result."""
        out = result.head(1).reset_index()
        for _, row in out.iterrows():
            click.echo(f"target: {row['target']}")
            click.echo(f"query: {row['query']}")
            click.echo(f"fident: {row['fident']}")
            click.echo(f"alnlen: {row['alnlen']}")
            click.echo(f"coverage: {row['coverage']}x")
            click.echo()

    def _filter_databases(self, db_name: str = None) -> dict:
        """Filter databases to search."""
        if db_name is None:
            return self.database_index

        if db_name not in self.database_index:
            raise ValueError(f"Database {db_name} not found in index!")

        return {db_name: self.database_index[db_name]}

    def run(
        self,
        reads: Path,
        n_reads: int = None,
        db_name: str = None,
        prefix: str = None,
    ) -> dict:
        """Run the complete search pipeline."""

        # Preprocess reads
        preprocessor = ReadPreprocessor(self.outdir / "tmp")
        subsampled_reads = preprocessor.preprocess_reads(reads, n_reads)

        # Create analyser
        analyser = ReadAnalyser(
            outdir=self.outdir,
            threads=self.threads,
            prefix=prefix or reads.stem,
            assembler=self.assembler,
        )

        # Determine which databases to search
        databases = self._filter_databases(db_name)

        # Run analysis
        results = {}
        for db_name, params in databases.items():
            result = analyser.analyse_short_reads(
                input_reads=subsampled_reads,
                marker=params.get("marker"),
                db=Path(params.get("db")),
                min_seq_id=params.get("min_seq_id"),
                min_aln_len=params.get("min_aln_len"),
            )
            if result is not None:
                results[db_name] = result

        # Display results
        for db_name, result in results.items():
            click.echo()
            click.echo(f"Top result for {db_name}:")
            self._print_result_summary(result)

        return results

    def cleanup(self) -> None:
        shutil.rmtree(self.outdir / "tmp")
