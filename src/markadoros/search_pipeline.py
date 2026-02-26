import shutil
from pathlib import Path

import click
import pandas as pd
from loguru import logger

from markadoros.assembler_runners import HifiasmRunner, SpadesRunner
from markadoros.contig_searcher import ContigSearcher
from markadoros.read_assembler import ReadAssembler
from markadoros.read_preprocessor import ReadPreprocessor


class SearchPipeline:
    """Orchestrates the barcode search workflow.

    Handles two input modes:
    1. Raw reads: preprocess -> assemble -> search contigs -> process results
    2. Pre-assembled contigs: search contigs -> process results
    """

    def __init__(
        self,
        outdir: Path,
        threads: int,
        database_index: dict,
        type: str,
        include_lineage: bool,
    ):
        self.outdir = Path(outdir)
        self.threads = threads
        self.database_index = database_index
        self.type = type
        self.include_lineage = include_lineage

        # Initialize assembler based on platform
        if type in ["sr", "short", "illumina"]:
            assembler = SpadesRunner(threads=self.threads, rna=False)
        elif type == "rnaseq":
            assembler = SpadesRunner(threads=self.threads, rna=True)
        elif type in ["pacbio_hifi", "pb"]:
            assembler = HifiasmRunner(threads=self.threads, ont=False)
        elif type in ["oxford_nanopore", "ont"]:
            assembler = HifiasmRunner(threads=self.threads, ont=True)
        elif type == "contigs":
            assembler = None
        else:
            raise ValueError(f"Unsupported input type: {type}")

        self.assembler = assembler

    def _print_result_summary(self, result: pd.DataFrame) -> None:
        """Print a summary of the top search result."""
        if result.empty:
            return

        out = result.head(1).reset_index()
        for _, row in out.iterrows():
            click.echo("")
            click.echo(f"contig: {row['target']}")
            click.echo(f"match_id: {row['seq_id']}")
            click.echo(f"match_taxon: {row['taxon']}")
            click.echo(f"fident: {row['fident']}")
            click.echo(f"alnlen: {row['alnlen']}")
            if "coverage" in out.columns:
                click.echo(f"coverage: {row['coverage']}x")

    def _filter_databases(self, db_name: str | None = None) -> dict:
        """Filter databases to search."""
        if db_name is None:
            return self.database_index

        if db_name not in self.database_index:
            raise ValueError(f"Database {db_name} not found in index!")

        return {db_name: self.database_index[db_name]}

    def _process_results(
        self,
        result: pd.DataFrame,
        marker: str,
        prefix: str,
        extract_coverage: bool = False,
        include_lineage: bool = False,
    ) -> pd.DataFrame:
        """Process and save search results.

        Extracts coverage from contig headers (if requested) and sorts by quality metrics.

        Args:
            result: DataFrame with search results
            marker: Marker name for output file
            prefix: Prefix for output file
            extract_coverage: Whether to extract coverage from target column
            include_lineage: Whether to extract lineage from target column
        """
        if result.empty:
            return result

        result["seq_id"] = result["query"].str.split("|").str[0]
        result["marker"] = result["query"].str.split("|").str[1]
        result["taxon"] = result["query"].str.split("|").str[2]

        if extract_coverage:
            result["coverage"] = (
                result["target"]
                .str.extract(r"cov_(\d+\.?\d*)", expand=False)
                .astype(float)
                .round(2)
            )
        else:
            result["coverage"] = pd.NA

        if include_lineage:
            result["lineage"] = result["query"].str.split("|").str[3]
        else:
            result["lineage"] = pd.NA

        if extract_coverage:
            result = result.sort_values(
                by=["coverage", "fident", "alnlen"], ascending=False
            )
        else:
            result = result.sort_values(by=["fident", "alnlen"], ascending=False)

        # Reorder columns: target, coverage, seq_id, marker, taxon, lineage, then rest
        desired_cols = ["target", "coverage", "seq_id", "marker", "taxon", "lineage"]

        # Add remaining columns in their original order
        remaining_cols = [
            col for col in result.columns if col not in desired_cols and col != "query"
        ]
        result = result[desired_cols + remaining_cols]

        result.to_csv(
            self.outdir / f"{prefix}.{marker}.result.tsv",
            index=False,
            sep="\t",
            na_rep="NA",
        )

        return result

    def run(
        self,
        input: Path,
        n_reads: int | None = None,
        db_name: str | None = None,
        prefix: str | None = None,
    ) -> dict:
        """Run the complete search pipeline.

        Args:
            reads: Path to reads file (FASTX or CRAM)
            n_reads: Optional number of reads to subsample
            db_name: Optional single database to search
            prefix: Output file prefix (defaults to reads stem)

        Returns:
            Dictionary mapping database names to result DataFrames
        """
        prefix = prefix or input.stem

        # Preprocess reads
        preprocessor = ReadPreprocessor(self.outdir / "tmp")
        subsampled_reads = preprocessor.preprocess_reads(input, n_reads)

        read_assembler = None
        if self.type != "contigs" and self.assembler is not None:
            # Assemble reads to contigs
            read_assembler = ReadAssembler(
                outdir=self.outdir,
                tmpdir=self.outdir / "tmp",
                threads=self.threads,
                prefix=prefix,
                assembler=self.assembler,
            )

        # Determine which databases to search
        databases = self._filter_databases(db_name)

        # Run analysis for each database
        results = {}
        for db_name, params in databases.items():
            if read_assembler is not None:
                contigs = read_assembler.assemble(
                    input_reads=subsampled_reads,
                    marker=params.get("marker"),
                    db=Path(params.get("db")),
                )
            else:
                contigs = input

            # Search contigs against database
            contig_searcher = ContigSearcher(
                tmpdir=self.outdir / "tmp",
                threads=self.threads,
            )

            if contigs is not None:
                result = contig_searcher.search_contigs(
                    contigs=contigs,
                    marker_db=Path(params.get("db")),
                    marker=params.get("marker"),
                    min_seq_id=params.get("min_seq_id"),
                    min_aln_len=params.get("min_aln_len"),
                )

            if result is None or result.empty:
                logger.warning(f"No results found for {db_name}!.")
                continue

            # Process results (extract coverage, sort, save)
            result = self._process_results(
                result,
                params.get("marker"),
                prefix,
                extract_coverage=(self.type != "contigs"),
                include_lineage=self.include_lineage,
            )
            results[db_name] = result

        # Display results
        for db_name, result in results.items():
            logger.info(f"Top result for {db_name}:")
            self._print_result_summary(result)

        return results

    def cleanup(self) -> None:
        """Remove temporary files."""
        shutil.rmtree(self.outdir / "tmp")
