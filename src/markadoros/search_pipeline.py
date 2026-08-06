from __future__ import annotations

import gzip
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
from markadoros.results_processor import ResultsProcessor

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
        synonyms: list[str],
        min_seq_id: float,
        min_aln_len: int,
        save_contigs: bool,
        fallback_reads: int,
    ):
        self.outdir = Path(outdir)
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)
        self.tmpdir = Path(tmpdir)
        self.threads = threads
        self.database_index = database_index
        self.db_to_tmpdir = db_to_tmpdir
        self.input_type = input_type
        self.expected_taxon = expected_taxon
        self.synonyms = synonyms
        self.all_taxon_names = [expected_taxon] + synonyms if expected_taxon else []
        self.min_seq_id = min_seq_id
        self.min_aln_len = min_aln_len
        self.save_contigs = save_contigs
        self.fallback_reads = fallback_reads

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

    def _setup_database(self, db: Path, prefix: str) -> Path:
        """Copy database to tmpdir to avoid overuse of disk."""
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
            tmpdir=self.tmpdir,
            threads=self.threads,
            prefix=prefix,
            assembler=self.assembler,
            input_type=self.input_type,
            fallback_reads=self.fallback_reads,
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
            Tuple of (n_aligned_reads, path_to_contigs)
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

    def _get_taxa_counts(self, taxa: list[str], taxon_json: Path) -> dict[str, int]:
        """Get count for specific taxa from the count JSON."""
        with gzip.open(taxon_json, "rt", encoding="utf-8") as f:
            data = json.load(f)
            return {taxon: data.get(taxon, 0) for taxon in taxa if taxon in data}

    def run(
        self,
        input: Path,
        db_name: str,
        n_reads: int | None = None,
        prefix: str | None = None,
    ) -> None:
        """Run the complete search pipeline.

        Orchestrates the entire workflow:
        1. Prepare workspace
        2. Setup preprocessing/assembly if needed
        3. Get or assemble contigs
        4. Search contigs against database
        5. Process and summarize results
        6. Save outputs

        Args:
            input: Path to input file (reads or contigs)
            db_name: Name of database to search
            n_reads: Optional number of reads to subsample
            prefix: Output file prefix (defaults to input stem)
        """
        # Prepare workspace
        self._prepare_workspace()
        output_prefix = prefix or input.stem

        # Validate database exists
        if db_name not in self.database_index:
            raise ValueError(f"Database {db_name} not found in index!")
        db_params = self.database_index[db_name]

        # Setup preprocessing and assembly if working with reads
        read_assembler, subsampled_reads = self._setup_preprocessing(
            input,
            n_reads,
            output_prefix,
        )

        # Setup database (copy to tmpdir if requested)
        if self.db_to_tmpdir:
            logger.info(f"Cloning database {db_name}...")
            db_path = self._setup_database(Path(db_params["db"]), output_prefix)
        else:
            db_path = Path(db_params["db"])

        # Get taxon count expectation
        taxa_counts = {}
        if self.expected_taxon:
            taxa_counts = self._get_taxa_counts(
                self.all_taxon_names, db_params["taxon_db"]
            )
            expected_taxon_count = taxa_counts.get(self.expected_taxon, 0)
            synonym_count = sum(taxa_counts.values()) - expected_taxon_count
            logger.info(
                f"There are {expected_taxon_count} possible records for {self.expected_taxon}, and {synonym_count} for synonyms."
            )

        marker = db_params["marker"]

        # Get contigs (either from input or assembly)
        n_aligned_reads, contigs = self._get_contigs(
            input=input,
            read_assembler=read_assembler,
            subsampled_reads=subsampled_reads,
            marker=marker if marker is not None else "",
            db_path=db_path,
        )

        if contigs is None:
            logger.warning(f"No contigs available for {db_name}!")
            return

        # Save contigs if requested
        if self.save_contigs:
            logger.info(
                f"Saving contigs to {self.outdir}/{output_prefix}.contigs.fasta"
            )
            shutil.copy(contigs, f"{self.outdir}/{output_prefix}.contigs.fasta")

        # Get assembly statistics
        contig_stats = self._get_contig_stats(contigs)
        logger.info(
            f"Assembly: {contig_stats['n']} contigs (longest {contig_stats['longest']}), with an N50 of {contig_stats['n50']}"
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
            return

        # Process results and generate summary
        results_processor = ResultsProcessor(
            result=result,
            contigs=contigs,
            input=input,
            marker=marker,
            database=str(db_params["db"]),
            expected_taxon=self.expected_taxon,
            synonyms=self.synonyms,
            input_type=self.input_type,
            n_reads=n_reads,
            n_aligned_reads=n_aligned_reads,
            contig_stats=contig_stats,
            assembler=self.assembler,
            taxa_counts=taxa_counts,
            extract_coverage=(self.input_type != InputType.CONTIGS),
            bins_taxa_path=db_params["bins_taxa"],
        )

        # Save summary
        results_processor.save_summary(self.outdir, output_prefix, marker)

        # Print top result summary
        logger.info(f"Top result for {db_name}:")
        results_processor.print_result_summary()

    def cleanup(self) -> None:
        """Remove temporary files."""
        shutil.rmtree(self.tmpdir)
