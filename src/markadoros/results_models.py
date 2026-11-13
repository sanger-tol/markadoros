import gzip
import json
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from pydantic import BaseModel

from markadoros.taxonomy import BinsTaxa


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

    def set_expected_taxon_fields(
        self,
        expected_taxon: str | None,
        synonyms: list[str] | None,
        bin_taxa: list[str] | None,
    ) -> None:
        """Classify a taxon against expected taxon and synonyms."""
        if self.taxon == expected_taxon:
            self.is_expected_taxon = True
            self.expected_taxon_match_type = "direct_match"
        elif synonyms and self.taxon in synonyms:
            self.is_expected_taxon = True
            self.expected_taxon_match_type = "synonym_match"
        elif bin_taxa and any(t == expected_taxon for t in bin_taxa):
            self.is_expected_taxon = True
            self.expected_taxon_match_type = "bold_bin_contains_expected_taxon"
        elif bin_taxa and synonyms and any(t in synonyms for t in bin_taxa):
            self.is_expected_taxon = True
            self.expected_taxon_match_type = "bold_bin_contains_synonym"
        else:
            self.is_expected_taxon = False
            self.expected_taxon_match_type = None

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
    taxon: str
    synonyms: list[str] | None
    counts: dict[str, int] | None


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
    marker: str | None
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

    def update_expected_species(
        self,
        expected_taxon: str,
        synonyms: list[str] | None,
        bin_taxa: BinsTaxa | None,
    ):
        for hit in self.iter_hit_records():
            hit.set_expected_taxon_fields(
                expected_taxon=expected_taxon,
                synonyms=synonyms,
                bin_taxa=bin_taxa.get_taxon(hit.target) if bin_taxa else None,
            )

        with gzip.open(self.input.database.parent / "taxon.json.gz", "rt") as fh:
            taxa = json.load(fh)

        self.input.expected_taxon = ExpectedTaxonInfo(
            taxon=expected_taxon,
            synonyms=synonyms,
            counts=taxa.get(expected_taxon, None),
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
