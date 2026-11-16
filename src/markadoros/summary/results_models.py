import gzip
import importlib.metadata
import json
import re
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
import pysam
from loguru import logger
from pydantic import BaseModel, field_validator

from markadoros.taxonomy.taxonomy import BinsTaxa
from markadoros.utils.input_types import InputType
from markadoros.utils.utils import extract_subsequence, get_mmseqs_version

if TYPE_CHECKING:
    from markadoros.assembly.assembler_runners import AssemblerRunner
    from markadoros.assembly.contigs import Contigs


class HitRecord(BaseModel):
    """One processed search hit (a row of the results table)."""

    target: str
    coverage: float | None = None
    seq_id: str
    marker: str
    taxon: str
    lineage: str
    is_expected_taxon: bool | None = None
    expected_taxon_match_type: str | None = None
    fident: float
    alnlen: int
    tstart: int
    tend: int
    evalue: float
    bits: float
    sequence: str

    def model_post_init(self, context: Any, /):
        if self.coverage is None:
            self.parse_coverage()

    def parse_coverage(self) -> None:
        match = re.search(r"cov_(\d+\.?\d*)", self.target)
        self.coverage = round(float(match.group(1)), 2) if match else None

    def to_tophit(self):
        """Return a TopHit object with a subset of the fields."""
        return TopHit(
            fident=self.fident,
            alnlen=self.alnlen,
            tstart=self.tstart,
            tend=self.tend,
            evalue=self.evalue,
            bits=self.bits,
            sequence=self.sequence,
        )


class TopHit(BaseModel):
    """Subset of a TopHit record for summary purposes."""

    fident: float
    alnlen: int
    tstart: int
    tend: int
    evalue: float
    bits: float
    sequence: str


class TaxonSummary(BaseModel):
    is_expected_taxon: bool | None
    expected_taxon_match_type: str | None
    n_hits: int
    fident_range: tuple[float, float]
    alnlen_range: tuple[int, int]
    top_hit: TopHit


class Summary(BaseModel):
    n_contigs_with_hits: int
    n_expected_taxon_hits: int | None
    n_synonym_hits: int | None
    top_result: HitRecord
    expected_taxon_top_result: HitRecord | None
    taxon_summary: dict[str, TaxonSummary]


class ExpectedTaxonInfo(BaseModel):
    """The taxon a run is searching for, its synonyms, and their counts in the database."""

    taxon: str
    synonyms: list[str] | None = None
    counts: dict[str, int] | None = None

    @field_validator("synonyms")
    @classmethod
    def normalise_synonyms(cls, synonyms: list[str] | None) -> list[str] | None:
        return synonyms or None

    @property
    def names(self) -> list[str]:
        """The expected taxon followed by its synonyms."""
        return [self.taxon, *(self.synonyms or [])]

    def annotate_hit(self, hit: HitRecord, bin_taxa: list[str] | None) -> None:
        """Set whether a hit's taxon matches the expected taxon or its synonyms."""
        synonyms = self.synonyms or []
        if hit.taxon == self.taxon:
            hit.is_expected_taxon = True
            hit.expected_taxon_match_type = "direct_match"
        elif synonyms and hit.taxon in synonyms:
            hit.is_expected_taxon = True
            hit.expected_taxon_match_type = "synonym_match"
        elif bin_taxa and any(t == self.taxon for t in bin_taxa):
            hit.is_expected_taxon = True
            hit.expected_taxon_match_type = "bold_bin_contains_expected_taxon"
        elif bin_taxa and synonyms and any(t in synonyms for t in bin_taxa):
            hit.is_expected_taxon = True
            hit.expected_taxon_match_type = "bold_bin_contains_synonym"
        else:
            hit.is_expected_taxon = False
            hit.expected_taxon_match_type = None


class ContigStats(BaseModel):
    """Assembly statistics for an assembly."""

    n: int
    size: int
    n50: int
    longest: int


class InputInfo(BaseModel):
    file: Path
    n_reads: int | None
    n_aligned_reads: int
    marker: str
    database: Path
    contig_stats: ContigStats | None
    expected_taxon: ExpectedTaxonInfo | None


class Tools(BaseModel):
    mmseqs: str | None
    assembler: dict[str, str | None] | None


class RunInfo(BaseModel):
    version: str
    tools: Tools


class ResultsReport(BaseModel):
    """Complete, validated output of a markadoros run."""

    input: InputInfo
    summary: Summary | None
    results: dict[str, list[HitRecord]] | None
    run_info: RunInfo

    @classmethod
    def from_search_results(
        cls,
        result: pd.DataFrame,
        contigs: "Contigs | None",
        input: Path,
        marker: str,
        database: str,
        expected_taxon: ExpectedTaxonInfo | None = None,
        input_type: InputType | None = None,
        n_reads: int | None = None,
        n_aligned_reads: int = 0,
        assembler: "AssemblerRunner | None" = None,
        extract_coverage: bool = False,
        bins_taxa_path: Path | None = None,
    ) -> "ResultsReport":
        """Process raw search results and build a validated report.

        Args:
            result: DataFrame with raw search results
            contigs: Contigs object for the contigs file, or None if no contigs were available
            input: Path to input file
            marker: Marker name
            database: Database path/name
            expected_taxon: Expected taxon with its synonyms and database counts (optional)
            input_type: Type of input (reads/contigs)
            n_reads: Number of reads processed
            n_aligned_reads: Number of reads with alignments
            assembler: Assembler instance for version info

            extract_coverage: Whether to extract coverage from contigs
            bins_taxa_path: Path to a gzipped JSON file containing the list of taxon names attached to each BOLD BIN

        Returns:
            ResultsReport with processed hits and summary
        """
        input_type = input_type or InputType.CONTIGS
        bins_taxa = BinsTaxa(bins_taxa_path)
        contig_stats = contigs._get_stats() if contigs is not None else None

        hits = cls._process_results(
            result,
            contigs.path if contigs is not None else None,
            extract_coverage,
            expected_taxon,
            bins_taxa,
        )
        summary, results = cls._summarise_results(hits, expected_taxon)

        return cls(
            input=InputInfo(
                file=input.resolve(),
                n_reads=n_reads if input_type != InputType.CONTIGS else None,
                n_aligned_reads=n_aligned_reads,
                marker=marker,
                database=database,
                contig_stats=contig_stats,
                expected_taxon=expected_taxon,
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

    @staticmethod
    def _process_results(
        result: pd.DataFrame,
        contigs: Path | None,
        extract_coverage: bool,
        expected_taxon: ExpectedTaxonInfo | None,
        bins_taxa: BinsTaxa,
    ) -> list[HitRecord]:
        """Process and enhance search results.

        Extracts marker info from headers, calculates coverage, and sorts results.

        Returns:
            List of HitRecords sorted by evalue
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
        result = result.sort_values(by=["evalue"], ascending=[False])
        records = result.to_dict(orient="records")

        results_list = [HitRecord(**record) for record in records]
        [record.parse_coverage() for record in results_list if extract_coverage]
        if expected_taxon:
            for record in results_list:
                expected_taxon.annotate_hit(record, bins_taxa.get_taxon(record.target))

        return results_list

    @staticmethod
    def _summarise_results(
        hits: list[HitRecord],
        expected_taxon: ExpectedTaxonInfo | None,
    ) -> tuple[Summary | None, dict[str, list[HitRecord]] | None]:
        """Summarise processed hits per taxon and per target contig.

        Returns:
            Tuple of (summary, results) for a non-empty hit list
        """
        if len(hits) == 0:
            return None, None

        all_taxon_names = expected_taxon.names if expected_taxon else []

        # Tally the number of hits per taxon
        found_taxa = {hit.taxon for hit in hits if hit.taxon is not None}

        # Top result - if expected taxon found, top result for taxon, otherwise overall top result
        top_result = hits[0]

        expected_taxon_top_result = None
        if expected_taxon and any(taxon in found_taxa for taxon in all_taxon_names):
            expected_taxon_top_result = next(
                hit for hit in hits if hit.taxon in all_taxon_names
            )

        # Summary for each taxon
        taxon_summary = {}
        for taxon in found_taxa:
            taxon_hits = [hit for hit in hits if hit.taxon == taxon]
            top_taxon_hit = taxon_hits[0]
            fidents = [hit.fident for hit in taxon_hits]
            alnlens = [hit.alnlen for hit in taxon_hits]
            taxon_summary[taxon] = TaxonSummary(
                is_expected_taxon=top_taxon_hit.is_expected_taxon,
                expected_taxon_match_type=top_taxon_hit.expected_taxon_match_type,
                n_hits=len(taxon_hits),
                fident_range=(min(fidents), max(fidents)),
                alnlen_range=(min(alnlens), max(alnlens)),
                top_hit=top_taxon_hit.to_tophit(),
            )

        # Get the number of hits for the expected taxon
        expected_taxon_counts_in_result = None
        synonym_taxon_counts_in_result = None
        if expected_taxon:
            synonyms = expected_taxon.synonyms or []
            expected_taxon_counts_in_result = sum(
                1 if hit.taxon == expected_taxon.taxon else 0 for hit in hits
            )
            synonym_taxon_counts_in_result = sum(
                1 if hit.taxon in synonyms else 0 for hit in hits
            )
            logger.info(
                f"Found {expected_taxon_counts_in_result} results for {expected_taxon.taxon} and {synonym_taxon_counts_in_result} for its synonyms."
            )

        summary = Summary(
            n_contigs_with_hits=len({hit.target for hit in hits}),
            n_expected_taxon_hits=expected_taxon_counts_in_result,
            n_synonym_hits=synonym_taxon_counts_in_result,
            top_result=top_result,
            expected_taxon_top_result=expected_taxon_top_result,
            taxon_summary=taxon_summary,
        )

        results = defaultdict(list)
        for item in hits:
            results[item.target].append(item)

        return summary, results

    def print_result_summary(self) -> None:
        """Print a summary of the top search result to logger."""
        if self.summary is None:
            return

        out = self.summary.top_result
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

    def update_expected_species(
        self,
        expected_taxon: ExpectedTaxonInfo,
        bin_taxa: BinsTaxa | None,
    ):
        for hit in self.iter_hit_records():
            expected_taxon.annotate_hit(
                hit, bin_taxa.get_taxon(hit.target) if bin_taxa else None
            )

        with gzip.open(self.input.database.parent / "taxon.json.gz", "rt") as fh:
            taxa = json.load(fh)

        self.input.expected_taxon = expected_taxon.model_copy(
            update={
                "counts": {
                    name: taxa[name] for name in expected_taxon.names if name in taxa
                }
            }
        )

        if self.results is not None and self.summary is not None:
            flat_records = sorted(
                (hit for hits in self.results.values() for hit in hits),
                key=lambda hit: hit.evalue,
            )
            expected_taxon_hits = [
                hit for hit in flat_records if hit.is_expected_taxon is True
            ]
            synonym_hits = [
                hit
                for hit in flat_records
                if hit.expected_taxon_match_type == "synonym_match"
            ]
            self.summary.n_expected_taxon_hits = len(expected_taxon_hits)
            self.summary.n_synonym_hits = len(synonym_hits)
            self.summary.expected_taxon_top_result = (
                expected_taxon_hits[0] if expected_taxon_hits else None
            )
            for taxon, taxon_summary in self.summary.taxon_summary.items():
                representative = next(hit for hit in flat_records if hit.taxon == taxon)
                taxon_summary.is_expected_taxon = representative.is_expected_taxon
                taxon_summary.expected_taxon_match_type = (
                    representative.expected_taxon_match_type
                )

    def iter_hit_records(self) -> Iterator[HitRecord]:
        """Yield every HitRecord in the summary and per-target results."""
        if self.summary is not None:
            yield self.summary.top_result
            if self.summary.expected_taxon_top_result is not None:
                yield self.summary.expected_taxon_top_result

        if self.results is not None:
            for hits in self.results.values():
                yield from hits
