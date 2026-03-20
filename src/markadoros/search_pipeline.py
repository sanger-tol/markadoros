from __future__ import annotations

import gzip
import importlib.metadata
import json
import shutil
from math import floor
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import pysam
from loguru import logger

from markadoros.assembler_runners import HifiasmRunner, SpadesRunner
from markadoros.contig_searcher import ContigSearcher
from markadoros.input_types import InputType
from markadoros.read_assembler import ReadAssembler
from markadoros.read_preprocessor import ReadPreprocessor
from markadoros.utils import extract_subsequence

if TYPE_CHECKING:
    from markadoros.assembler_runners import AssemblerRunner


class SearchPipeline:
    """Orchestrates the barcode search workflow.

    Handles two input modes:
    1. Raw reads: preprocess -> assemble -> search contigs -> process results
    2. Pre-assembled contigs: search contigs -> process results
    """

    def __init__(
        self,
        outdir: Path,
        tmpdir: Path,
        threads: int,
        database_index: dict[str, dict],
        db_to_tmpdir: bool,
        input_type: InputType,
        expected_taxon: str | None,
        min_seq_id: float,
        min_aln_len: int,
    ):
        self.outdir = Path(outdir)
        self.tmpdir = Path(tmpdir)
        self.threads = threads
        self.database_index = database_index
        self.db_to_tmpdir = db_to_tmpdir
        self.input_type = input_type
        self.expected_taxon = expected_taxon
        self.min_seq_id = min_seq_id
        self.min_aln_len = min_aln_len

        # Initialize assembler based on platform
        assembler: AssemblerRunner | None
        if input_type == InputType.SHORT_READ:
            assembler = SpadesRunner(threads=self.threads, rna=False)
        elif input_type == InputType.RNASEQ:
            assembler = SpadesRunner(threads=self.threads, rna=True)
        elif input_type == InputType.PACBIO:
            assembler = HifiasmRunner(threads=self.threads, ont=False)
        elif input_type == InputType.ONT:
            assembler = HifiasmRunner(threads=self.threads, ont=True)
        elif input_type == InputType.CONTIGS:
            assembler = None
        else:
            raise ValueError(f"Unsupported input type: {input_type}")

        self.assembler = assembler

    def _prepare_workspace(self) -> None:
        """Prepare the workspace, cleaning up any leftover tmpdir from previous runs."""
        if self.tmpdir.exists():
            logger.warning(
                f"Cleaning up existing tmpdir from previous run: {self.tmpdir}"
            )
            shutil.rmtree(self.tmpdir)
            self.tmpdir.mkdir(parents=True, exist_ok=True)
        else:
            self.tmpdir.mkdir(parents=True, exist_ok=True)

    def _print_result_summary(self, result: pd.DataFrame) -> None:
        """Print a summary of the top search result."""
        if result.empty:
            return

        if self.expected_taxon and self.expected_taxon in result["taxon"].values:
            out = result[result["taxon"] == self.expected_taxon].head(1).reset_index()
        else:
            out = result.head(1).reset_index()
        for _, row in out.iterrows():
            logger.info(f"""\n
            \tcontig: {row["target"]}
            \tmatch_id: {row["seq_id"]}
            \tmatch_taxon: {row["taxon"]}
            \tfident: {row["fident"]}
            \talnlen: {row["alnlen"]}
            \tcoverage: {row["coverage"]}x
            """)

    def _filter_databases(self, db_name: str) -> tuple[str, dict]:
        """Filter databases to search."""
        if db_name not in self.database_index:
            raise ValueError(f"Database {db_name} not found in index!")

        return db_name, self.database_index[db_name]

    def _process_results(
        self,
        result: pd.DataFrame,
        contigs: Path | None,
        extract_coverage: bool = False,
    ) -> pd.DataFrame:
        """Process and save search results.

        Extracts coverage from contig headers (if requested) and sorts by quality metrics.

        Args:
            result: DataFrame with search results
            extract_coverage: Whether to extract coverage from target column
            include_lineage: Whether to extract lineage from target column
        """
        if result.empty:
            return result

        if contigs is None:
            raise FileNotFoundError("Error: Contigs file not found!")

        sequences = {}
        with pysam.FastxFile(str(contigs), persist=False) as fh:
            for record in fh:
                if record.sequence is not None:
                    sequences[record.name] = record.sequence.upper()

        result["seq_id"] = result["query"].str.split("|").str[0]
        result["marker"] = result["query"].str.split("|").str[1]
        result["taxon"] = result["query"].str.split("|").str[2].str.replace("_", " ")
        result["lineage"] = result["query"].str.split("|").str[3]
        result["sequence"] = result.apply(
            lambda row: extract_subsequence(row, sequences), axis=1
        )

        if extract_coverage:
            result["coverage"] = (
                result["target"]
                .str.extract(r"cov_(\d+\.?\d*)", expand=False)
                .astype(float)
                .round(2)
            )
        else:
            result["coverage"] = pd.NA

        if extract_coverage:
            result = result.sort_values(
                by=["coverage", "fident", "alnlen"], ascending=False
            )
        else:
            result = result.sort_values(by=["target", "bits"], ascending=[True, False])

        # Reorder columns: target, coverage, seq_id, marker, taxon, lineage, then rest
        desired_cols = ["target", "coverage", "seq_id", "marker", "taxon", "lineage"]

        # Add remaining columns in their original order
        remaining_cols = [
            col for col in result.columns if col not in desired_cols and col != "query"
        ]
        result = result[desired_cols + remaining_cols]

        return result

    def _summarise_result(
        self,
        input: Path,
        result: pd.DataFrame,
        taxon_count: int | None,
        marker: str | None,
        database: str,
        n_reads: int | None,
        n_aligned_reads: int,
        contig_stats: dict | None = None,
    ) -> dict:
        """
        Summarise a results dataframe, and if provided an expected taxon, check if
        it was found and how many results were found.
        """
        ## Tally the number of hits per taxon
        found_taxon_counts = result["taxon"].value_counts()

        ## Top result - if expected taxon found, top result for taxon, otherwise overall top result
        if self.expected_taxon and self.expected_taxon in found_taxon_counts:
            top_result = (
                result[result["taxon"] == self.expected_taxon]
                .iloc[0]
                .to_dict(orient="records")
            )
        else:
            top_result = result.iloc[0].to_dict(orient="records")

        ## Summary for each taxon
        taxon_summary = {}
        for taxon in found_taxon_counts.index:
            taxon_hits = result[result["taxon"] == taxon]
            top_taxon_hit = taxon_hits.iloc[0]
            taxon_summary[taxon] = {
                "expected_taxon": True if self.expected_taxon == taxon else False,
                "n_hits": len(taxon_hits),
                "fident_range": [
                    float(taxon_hits["fident"].min()),
                    float(taxon_hits["fident"].max()),
                ],
                "alnlen_range": [
                    int(taxon_hits["alnlen"].min()),
                    int(taxon_hits["alnlen"].max()),
                ],
                "top_hit": {
                    "fident": float(top_taxon_hit["fident"]),
                    "alnlen": int(top_taxon_hit["alnlen"]),
                    "tstart": int(top_taxon_hit["tstart"]),
                    "tend": int(top_taxon_hit["tend"]),
                    "evalue": float(top_taxon_hit["evalue"]),
                    "bits": int(top_taxon_hit["bits"]),
                    "sequence": str(top_taxon_hit["sequence"]),
                },
            }

        ## Get the number of hits for the expected taxon
        expected_taxon_counts_in_result = None
        if self.expected_taxon:
            expected_taxon_counts_in_result = int(
                found_taxon_counts.get(self.expected_taxon, 0)
            )
            logger.info(
                f"Found {expected_taxon_counts_in_result} results for {self.expected_taxon}!"
            )

        expectation = (
            {
                "taxon": self.expected_taxon,
                "counts": {"available_sequences": taxon_count},
            }
            if self.expected_taxon
            else {}
        )

        summary = {
            "n_contigs_with_hits": int(result["target"].nunique()),
            "n_expected_taxon_hits": expected_taxon_counts_in_result,
            "top_result": top_result,
            "taxon_summary": taxon_summary,
        }

        results = {
            target: group.to_dict(orient="records")
            for target, group in result.groupby("target")
        }

        output = {
            "input": {
                "file": str(input.resolve()),
                "n_reads": n_reads if self.input_type != InputType.CONTIGS else None,
                "n_aligned_reads": n_aligned_reads,
                "marker": marker,
                "database": database,
                "contig_stats": contig_stats,
                "expected_taxon": expectation,
            },
            "summary": summary,
            "results": results,
            "run_info": {
                "version": importlib.metadata.version("markadoros"),
            },
        }

        return output

    def _get_taxon_count(self, taxon_json: Path) -> int:
        """Get count for a specific species from the count JSON."""
        with gzip.open(taxon_json, "rt", encoding="utf-8") as f:
            data = json.load(f)
            return data.get(self.expected_taxon, 0)

    def _setup_database(self, db: Path, prefix: str) -> Path:
        """
        Copy database to tmpdir to avoid overuse of disk
        """
        out_db_path = self.tmpdir / "marker_db" / db.parent.name
        out_db_path.mkdir(parents=True, exist_ok=True)

        for file in db.parent.glob("*"):
            if not file.is_dir():
                shutil.copy2(file, out_db_path / file.name)

        return (out_db_path / db.name).resolve()

    def _get_contig_stats(self, contigs: Path) -> dict[str, int]:
        """Calculate assembly statistics from contigs.

        Args:
            contigs: Path to contigs FASTA file

        Returns:
            Dictionary with assembly statistics (n, size, n50, longest)
        """
        count = 0
        lengths: list[int] = []
        with pysam.FastxFile(str(contigs)) as asm:
            for record in asm:
                count += 1
                if record.sequence is not None:
                    lengths.append(len(record.sequence))

        if not lengths:
            logger.warning("No contigs found in assembly file!")
            return {"n": 0, "size": 0, "n50": 0, "longest": 0}

        asm_size = sum(lengths)
        longest = max(lengths)

        lengths.sort(reverse=True)
        cumulative_sum = 0
        n50 = 0
        for length in lengths:
            cumulative_sum += length
            if cumulative_sum >= floor(asm_size / 2):
                n50 = length
                break

        return {"n": count, "size": asm_size, "n50": n50, "longest": longest}

    def _setup_preprocessing(
        self, input: Path, n_reads: int | None, prefix: str
    ) -> tuple[ReadAssembler | None, Path | None]:
        """Setup read preprocessing and assembly pipeline if needed.

        Returns:
            Tuple of (read_assembler, subsampled_reads) or (None, None) if not applicable
        """
        if self.input_type == InputType.CONTIGS or self.assembler is None:
            return None, None

        preprocessor = ReadPreprocessor(self.tmpdir)
        subsampled_reads = preprocessor.preprocess_reads(input, n_reads)

        read_assembler = ReadAssembler(
            outdir=self.outdir,
            tmpdir=self.tmpdir,
            threads=self.threads,
            prefix=prefix,
            assembler=self.assembler,
            input_type=self.input_type,
        )

        return read_assembler, subsampled_reads

    def _get_contigs(
        self,
        input: Path,
        read_assembler: ReadAssembler | None,
        subsampled_reads: Path | None,
        marker: str,
        db_path: Path,
    ) -> tuple[int, Path | None]:
        """Get contigs either from assembly or input.

        Returns:
            Path to contigs file, or None if assembly failed
        """
        if read_assembler is not None and subsampled_reads is not None:
            return read_assembler.assemble(
                input_reads=subsampled_reads,
                marker=marker,
                db=db_path,
            )
        else:
            return 0, input

    def _search_contigs(
        self,
        contigs: Path,
        db_path: Path,
        marker: str,
        min_seq_id: float | None,
        min_aln_len: int | None,
    ) -> pd.DataFrame | None:
        """Search contigs against the marker database.

        Returns:
            DataFrame with search results, or None if search failed
        """
        contig_searcher = ContigSearcher(
            tmpdir=self.tmpdir,
            threads=self.threads,
        )

        return contig_searcher.search_contigs(
            contigs=contigs,
            marker_db=db_path,
            marker=marker,
            min_seq_id=min_seq_id if min_seq_id is not None else 0.0,
            min_aln_len=min_aln_len if min_aln_len is not None else 0,
        )

    def _save_summary(
        self,
        summary: dict,
        prefix: str,
        marker: str | None,
    ) -> None:
        """Save summary results to JSON file."""
        output_path = self.outdir / f"{prefix}.{marker}.summary.json"
        with open(output_path, "w") as f:
            json.dump(summary, f, indent=4)

    def _process_database(
        self,
        db_name: str,
        params: dict,
        input: Path,
        read_assembler: ReadAssembler | None,
        subsampled_reads: Path | None,
        prefix: str,
        n_reads: int | None,
    ) -> pd.DataFrame | None:
        """Process a single database: assemble (if needed), search, and save results.

        Returns:
            Processed results DataFrame, or None if no results found
        """
        if self.db_to_tmpdir:
            logger.info(f"Cloning database {db_name}...")
            db_path = self._setup_database(Path(params["db"]), prefix)
        else:
            db_path = Path(params["db"])

        # Get taxon count expectation
        taxon_count = None
        if self.expected_taxon:
            taxon_count = self._get_taxon_count(params["taxon_db"])
            logger.info(
                f"There are {taxon_count} possible records for {self.expected_taxon}!"
            )

        marker = params["marker"]

        n_aligned_reads, contigs = self._get_contigs(
            input=input,
            read_assembler=read_assembler,
            subsampled_reads=subsampled_reads,
            marker=marker if marker is not None else "",
            db_path=db_path,
        )

        if contigs is None:
            logger.warning(f"No contigs available for {db_name}!")
            return None

        contigs_stats = self._get_contig_stats(contigs)
        logger.info(
            f"Assembly: {contigs_stats['n']} contigs (longest {contigs_stats['longest']}), with an N50 of {contigs_stats['n50']}"
        )

        # Search contigs against database
        result = self._search_contigs(
            contigs=contigs,
            db_path=db_path,
            marker=marker if marker is not None else "",
            min_seq_id=self.min_seq_id,
            min_aln_len=self.min_aln_len,
        )

        # Handle empty results
        if result is None or result.empty:
            logger.warning(f"No results found for {db_name}!")
            return None

        # Process results (extract coverage, sort, etc.)
        result = self._process_results(
            result=result,
            extract_coverage=(self.input_type != InputType.CONTIGS),
            contigs=contigs,
        )

        # Generate and save summary
        summary = self._summarise_result(
            input=input,
            result=result,
            taxon_count=taxon_count,
            marker=marker,
            database=str(params["db"]),
            n_reads=n_reads,
            n_aligned_reads=n_aligned_reads,
            contig_stats=contigs_stats,
        )
        self._save_summary(summary, prefix, marker)

        return result

    def run(
        self,
        input: Path,
        db_name: str,
        n_reads: int | None = None,
        prefix: str | None = None,
    ) -> None:
        """Run the complete search pipeline.

        Args:
            input: Path to input file (reads or contigs)
            n_reads: Optional number of reads to subsample
            db_name: Optional single database to search
            prefix: Output file prefix (defaults to input stem)

        Returns:
            Dictionary mapping database names to result DataFrames
        """
        # Prepare workspace (clean up any leftover tmpdir)
        self._prepare_workspace()

        output_prefix = prefix or input.stem

        # Setup preprocessing and assembly if working with reads
        read_assembler, subsampled_reads = self._setup_preprocessing(
            input, n_reads, output_prefix
        )

        # Extract relevant DB parameters
        db_name, params = self._filter_databases(db_name)

        result = self._process_database(
            db_name=db_name,
            params=params,
            input=input,
            read_assembler=read_assembler,
            subsampled_reads=subsampled_reads,
            prefix=output_prefix,
            n_reads=n_reads,
        )

        if result is None:
            return None

        logger.info(f"Top result for {db_name}:")
        self._print_result_summary(result)

        return None

    def cleanup(self) -> None:
        """Remove temporary files."""
        shutil.rmtree(self.tmpdir)
