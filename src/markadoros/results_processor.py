import importlib.metadata
import json
from pathlib import Path
from typing import TYPE_CHECKING
import gzip

import pandas as pd
import pysam
from loguru import logger

from markadoros.input_types import InputType
from markadoros.utils import extract_subsequence, get_mmseqs_version

if TYPE_CHECKING:
    from markadoros.assembler_runners import AssemblerRunner


class ResultsProcessor:
    """Processes search results and generates summaries.

    Automatically processes results on initialization and stores them
    along with a generated summary for later access.
    """

    def __init__(
        self,
        result: pd.DataFrame,
        contigs: Path,
        input: Path,
        marker: str | None,
        database: str,
        expected_taxon: str | None = None,
        synonyms: list[str] | None = None,
        input_type: InputType | None = None,
        n_reads: int | None = None,
        n_aligned_reads: int = 0,
        contig_stats: dict | None = None,
        assembler: "AssemblerRunner | None" = None,
        taxa_counts: dict[str, int] | None = None,
        extract_coverage: bool = False,
        bins_taxa_path: Path | None = None
    ):
        """Initialize processor with results and parameters.

        Automatically processes results and generates summary on initialization.

        Args:
            result: DataFrame with raw search results
            contigs: Path to contigs file
            input: Path to input file
            marker: Marker name
            database: Database path/name
            expected_taxon: Expected taxon name (optional)
            synonyms: List of synonym taxon names (optional)
            input_type: Type of input (reads/contigs)
            n_reads: Number of reads processed
            n_aligned_reads: Number of reads with alignments
            contig_stats: Assembly statistics
            assembler: Assembler instance
            taxa_counts: Expected taxon counts from database
            extract_coverage: Whether to extract coverage from contigs
            bins_taxa_path: Path to a gzipped JSON file containing the list of taxon names attached to each BOLD BIN
        """
        self.expected_taxon = expected_taxon
        self.synonyms = synonyms or []
        self.all_taxon_names = (
            [expected_taxon] + self.synonyms if expected_taxon else []
        )
        self.input_type = input_type or InputType.CONTIGS

        self.bins_taxa = self._load_bins_taxa(bins_taxa_path)

        # Process results
        self.result = self._process_results(result, contigs, extract_coverage)

        # Generate summary
        self.summary = self._generate_summary(
            input=input,
            taxa_counts=taxa_counts,
            marker=marker,
            database=database,
            n_reads=n_reads,
            n_aligned_reads=n_aligned_reads,
            contig_stats=contig_stats,
            assembler=assembler,
        )

    def _load_bins_taxa(self, path) -> dict[str, list[str]]:
        if not path:
            return {}
        with gzip.open(path, "rt") as fh:
            return json.load(fh)

    def _get_bin_taxa(self, seq_id: str) -> list[str]:
        parts = seq_id.split("/")
        if len(parts) >= 2:
            bin_id = parts[1]
            return self.bins_taxa.get(bin_id, [])
        return []

    def _get_expected_taxon_match_type(self, taxon, bin_taxa):
        """Classify a taxon against expected taxon and synonyms."""
        if taxon == self.expected_taxon:
            return "direct_match"
        if taxon in self.synonyms:
            return "synonym_match"
        if bin_taxa and any(t == self.expected_taxon for t in bin_taxa):
            return "bold_bin_contains_expected_taxon"
        if bin_taxa and any(t in self.synonyms for t in bin_taxa):
            return "bold_bin_contains_synonym"
        return "none"

    def _is_expected_taxon(self, taxon, bin_taxa):
        """Classify a taxon against expected taxon and synonyms."""
        if taxon == self.expected_taxon:
            return "true"
        if taxon in self.synonyms:
            return "true"
        if bin_taxa and any(t == self.expected_taxon for t in bin_taxa):
            return "true"
        if bin_taxa and any(t in self.synonyms for t in bin_taxa):
            return "true"
        return "false"

    def _process_results(
        self,
        result: pd.DataFrame,
        contigs: Path,
        extract_coverage: bool = False,
    ) -> pd.DataFrame:
        """Process and enhance search results.

        Extracts marker info from headers, calculates coverage, and sorts results.

        Args:
            result: DataFrame with search results
            contigs: Path to contigs file
            extract_coverage: Whether to extract coverage from target column

        Returns:
            Processed results DataFrame
        """
        if result.empty:
            return result

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
        result["bin_taxa"] = result["seq_id"].apply(
            lambda seq_id: self._get_bin_taxa(seq_id)
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

        result = result.sort_values(
            by=["bits"], ascending=[False]
        )

        if self.expected_taxon:
            result["is_expected_taxon"] = result.apply(
                lambda row: self._is_expected_taxon(row["taxon"], row["bin_taxa"]),
                axis=1
            )
            result["expected_taxon_match_type"] = result.apply(
                lambda row: self._get_expected_taxon_match_type(row["taxon"], row["bin_taxa"]),
                axis=1
            )
        result = result.drop(columns=["bin_taxa"])

        # Reorder columns: target, coverage, seq_id, marker, taxon, lineage, then rest
        desired_cols = [
            "target",
            "coverage",
            "seq_id",
            "marker",
            "taxon",
            "lineage",
        ]

        # Only include is_expected_taxon if it was created
        if "is_expected_taxon" in result.columns:
            desired_cols.append("is_expected_taxon")
            desired_cols.append("expected_taxon_match_type")

        # Add remaining columns in their original order
        remaining_cols = [
            col for col in result.columns if col not in desired_cols and col != "query"
        ]
        result = result[desired_cols + remaining_cols]

        return result

    def _generate_summary(
        self,
        input: Path,
        taxa_counts: dict[str, int] | None,
        marker: str | None,
        database: str,
        n_reads: int | None,
        n_aligned_reads: int,
        contig_stats: dict | None = None,
        assembler: "AssemblerRunner | None" = None,
    ) -> dict:
        """Generate summary for processed results.

        Args:
            input: Input file path
            taxa_counts: Expected taxon counts
            marker: Marker name
            database: Database path
            n_reads: Number of reads processed
            n_aligned_reads: Number of aligned reads
            contig_stats: Assembly statistics
            assembler: Assembler instance for version info

        Returns:
            Dictionary with complete summary
        """
        if self.result.empty:
            return {}

        # Statistics regarding the expected taxon, if provided
        expectation = (
            {
                "taxon": self.expected_taxon,
                "synonyms": self.synonyms if self.synonyms else None,
                "counts": taxa_counts,
            }
            if self.expected_taxon
            else {}
        )

        # Tally the number of hits per taxon
        found_taxon_counts = self.result["taxon"].value_counts()

        # Top result - if expected taxon found, top result for taxon, otherwise overall top result
        top_result = self.result.iloc[0].to_dict()

        expected_taxon_top_result = {}
        if self.expected_taxon and any(
            taxon in found_taxon_counts.index for taxon in self.all_taxon_names
        ):
            expected_taxon_top_result = (
                self.result[self.result["taxon"].isin(self.all_taxon_names)]
                .iloc[0]
                .to_dict()
            )

        # Summary for each taxon
        taxon_summary = {}
        for taxon in found_taxon_counts.index:
            taxon_hits = self.result[self.result["taxon"] == taxon]
            top_taxon_hit = taxon_hits.iloc[0]
            taxon_summary[taxon] = {
                "is_expected_taxon": self._is_expected_taxon(taxon, None),
                "expected_taxon_match_type": self._get_expected_taxon_match_type(taxon, None),
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

        # Get the number of hits for the expected taxon
        expected_taxon_counts_in_result = None
        synonym_taxon_counts_in_result = None
        if self.expected_taxon:
            expected_taxon_counts_in_result = int(
                found_taxon_counts.get(self.expected_taxon, 0)
            )
            synonym_taxon_counts_in_result = int(
                sum([found_taxon_counts.get(taxon, 0) for taxon in self.synonyms])
            )
            logger.info(
                f"Found {expected_taxon_counts_in_result} results for {self.expected_taxon} and {synonym_taxon_counts_in_result} for its synonyms."
            )

        summary = {
            "n_contigs_with_hits": int(self.result["target"].nunique()),
            "n_expected_taxon_hits": expected_taxon_counts_in_result,
            "n_synonym_hits": synonym_taxon_counts_in_result,
            "top_result": top_result,
            "expected_taxon_top_result": expected_taxon_top_result,
            "taxon_summary": taxon_summary,
        }

        results = {
            target: group.to_dict(orient="records")
            for target, group in self.result.groupby("target")
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
                "tools": {
                    "mmseqs": get_mmseqs_version(),
                    "assembler": {assembler.name: assembler._get_version()}
                    if assembler is not None
                    else None,
                },
            },
        }

        return output

    def get_summary(self) -> dict:
        """Get the generated summary.

        Returns:
            Dictionary with complete summary information
        """
        return self.summary

    def print_result_summary(self) -> None:
        """Print a summary of the top search result to logger."""
        if self.result.empty:
            return

        out = self.result.head(1).reset_index()

        for _, row in out.iterrows():
            is_expected_taxon = row.get("is_expected_taxon", "unknown")
            logger.info(
                f"""\n
            \tcontig: {row["target"]}
            \tmatch_id: {row["seq_id"]}
            \tmatch_taxon: {row["taxon"]}
            \tfident: {row["fident"]}
            \talnlen: {row["alnlen"]}
            \tcoverage: {row["coverage"]}x
            \tis_expected_taxon: {is_expected_taxon}
            """
            )

    def save_summary(self, outdir: Path, prefix: str, marker: str | None) -> None:
        """Save summary to JSON file.

        Args:
            outdir: Output directory
            prefix: Output file prefix
            marker: Marker name (included in filename)
        """
        output_path = outdir / f"{prefix}.{marker}.summary.json"
        with open(output_path, "w") as f:
            json.dump(self.summary, f, indent=4)
        logger.info(f"Saved summary to {output_path}")
