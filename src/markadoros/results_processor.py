import importlib.metadata
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import pysam
from loguru import logger

from markadoros.input_types import InputType
from markadoros.results_models import (
    ContigStats,
    ExpectedTaxonInfo,
    HitRecord,
    InputInfo,
    ResultsReport,
    RunInfo,
    Summary,
    TaxonSummary,
    Tools,
)
from markadoros.taxonomy import BinsTaxa
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
        contigs: Path | None,
        input: Path,
        marker: str | None,
        database: str,
        expected_taxon: str | None = None,
        synonyms: list[str] | None = None,
        input_type: InputType | None = None,
        n_reads: int | None = None,
        n_aligned_reads: int = 0,
        contig_stats: ContigStats | None = None,
        assembler: "AssemblerRunner | None" = None,
        taxa_counts: dict[str, int] | None = None,
        extract_coverage: bool = False,
        bins_taxa_path: Path | None = None,
    ):
        """Initialize processor with results and parameters.

        Automatically processes results and generates summary on initialization.

        Args:
            result: DataFrame with raw search results
            contigs: Path to contigs file, or None if no contigs were available
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

        self.bins_taxa = BinsTaxa(bins_taxa_path)

        # Process results
        self.result = self._process_results(result, contigs, extract_coverage)

        # Generate summary
        self.report = self._generate_summary(
            input=input,
            taxa_counts=taxa_counts,
            marker=marker,
            database=database,
            n_reads=n_reads,
            n_aligned_reads=n_aligned_reads,
            contig_stats=contig_stats,
            assembler=assembler,
        )

    def _process_results(
        self,
        result: pd.DataFrame,
        contigs: Path | None,
        extract_coverage: bool = False,
    ) -> list[HitRecord]:
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
            return []

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
        records = result.to_dict(orient="records")

        results_list = [HitRecord(**record) for record in records]
        [record.parse_coverage() for record in results_list if extract_coverage]
        if self.expected_taxon:
            for record in results_list:
                record.set_expected_taxon_fields(
                    expected_taxon=self.expected_taxon,
                    synonyms=self.synonyms,
                    bin_taxa=self.bins_taxa.get_taxon(record.target),
                )

        return results_list

    def _generate_summary(
        self,
        input: Path,
        taxa_counts: dict[str, int] | None,
        marker: str | None,
        database: str,
        n_reads: int | None,
        n_aligned_reads: int,
        contig_stats: ContigStats | None = None,
        assembler: "AssemblerRunner | None" = None,
    ) -> ResultsReport:
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
            Validated ResultsReport with complete summary
        """
        # Statistics regarding the expected taxon, if provided
        expectation = (
            ExpectedTaxonInfo(
                taxon=self.expected_taxon,
                synonyms=self.synonyms if self.synonyms else None,
                counts=taxa_counts,
            )
            if self.expected_taxon
            else None
        )

        summary, results = self._summarise_results()

        report = ResultsReport(
            input=InputInfo(
                file=str(input.resolve()),
                n_reads=n_reads if self.input_type != InputType.CONTIGS else None,
                n_aligned_reads=n_aligned_reads,
                marker=marker,
                database=database,
                contig_stats=contig_stats,
                expected_taxon=expectation or None,
            ),
            summary=summary,
            results=results,
            run_info=RunInfo(
                version=importlib.metadata.version("markadoros"),
                tools=Tools(
                    mmseqs=get_mmseqs_version(),
                    assembler={assembler.name: assembler._get_version()}
                    if assembler is not None
                    else None,
                ),
            ),
        )

        return report

    def _summarise_results(
        self,
    ) -> tuple[Summary | None, dict[str, list[HitRecord]] | None]:
        """Summarise processed results per taxon and per target contig.

        Returns:
            Tuple of (summary, results) for a non-empty result set
        """
        if len(self.result) == 0:
            return None, None

        # Tally the number of hits per taxon
        found_taxa = {hit.taxon for hit in self.result if hit.taxon is not None}

        # Top result - if expected taxon found, top result for taxon, otherwise overall top result
        top_result = self.result[0]

        expected_taxon_top_result = None
        if self.expected_taxon and any(
            taxon in found_taxa for taxon in self.all_taxon_names
        ):
            expected_taxon_top_result = next(
                hit for hit in self.result if hit.taxon in self.all_taxon_names
            )

        # Summary for each taxon
        taxon_summary = {}
        for taxon in found_taxa:
            taxon_hits = [hit for hit in self.result if hit.taxon == taxon]
            top_taxon_hit = taxon_hits[0]
            fidents = [hit.fident for hit in taxon_hits]
            alnlens = [hit.alnlen for hit in taxon_hits]
            taxon_summary[taxon] = TaxonSummary(
                is_expected_taxon=top_taxon_hit.is_expected_taxon
                if top_taxon_hit
                else None,
                expected_taxon_match_type=top_taxon_hit.expected_taxon_match_type
                if top_taxon_hit
                else None,
                n_hits=len(taxon_hits),
                fident_range=(min(fidents), max(fidents)),
                alnlen_range=(min(alnlens), max(alnlens)),
                top_hit=top_taxon_hit.to_tophit(),
            )

        # Get the number of hits for the expected taxon
        expected_taxon_counts_in_result = None
        synonym_taxon_counts_in_result = None
        if self.expected_taxon:
            expected_taxon_counts_in_result = sum(
                1 if hit.taxon == self.expected_taxon else 0 for hit in self.result
            )
            synonym_taxon_counts_in_result = sum(
                1 if hit.taxon in self.synonyms else 0 for hit in self.result
            )
            logger.info(
                f"Found {expected_taxon_counts_in_result} results for {self.expected_taxon} and {synonym_taxon_counts_in_result} for its synonyms."
            )

        summary = Summary(
            n_contigs_with_hits=len({hit.target for hit in self.result}),
            n_expected_taxon_hits=expected_taxon_counts_in_result,
            n_synonym_hits=synonym_taxon_counts_in_result,
            top_result=top_result,
            expected_taxon_top_result=expected_taxon_top_result,
            taxon_summary=taxon_summary,
        )

        results = defaultdict(list)
        for item in self.result:
            results[item.target].append(item)

        return summary, results

    def get_summary(self) -> ResultsReport:
        """Get the generated summary.

        Returns:
            ResultsReport with complete summary information
        """
        return self.report

    def print_result_summary(self) -> None:
        """Print a summary of the top search result to logger."""
        if len(self.result) == 0:
            return

        if out := self.result[0] if self.result else None:
            logger.info(
                f"""\n
            \tcontig: {out.target}
            \tmatch_id: {out.seq_id}
            \tmatch_taxon: {out.taxon}
            \tfident: {out.fident}
            \talnlen: {out.alnlen}
            \tcoverage: {out.coverage}x
            \tis_expected_taxon: {out.is_expected_taxon}
            \texpected_taxon_match_type: {out.expected_taxon_match_type}
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
        output_path.write_text(self.report.model_dump_json(indent=4))
        logger.info(f"Saved summary to {output_path}")
