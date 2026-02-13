import gzip
import shutil
from pathlib import Path

import pysam
from pymmseqs.config import CreateDBConfig as CreateMMSeqsDBConfig
from pymmseqs.config import EasySearchConfig as MMSeqsEasySearchConfig
from pymmseqs.parsers import EasySearchParser as MMSeqsEasySearchParser

from markadoros.spades import SpadesRunner


class ShortReadAnalyser:
    def __init__(
        self,
        outdir: Path,
        rna: bool,
        threads: int,
    ):
        self.outdir = outdir
        self.rna = rna
        self.threads = threads
        self.marker = None

        self.aligned_reads = None
        self.assembly_db = None

    def _get_marker_from_db(self, db: Path) -> str:
        """Extract marker name from database path."""
        return db.parent.name

    def _extract_reads(
        self, search_result: MMSeqsEasySearchParser, input_reads: Path
    ) -> Path:
        """Extract aligned reads to file. Returns path to output file."""
        aligned_reads = {read["query"] for read in search_result.to_gen()}
        outfile = self.outdir / self.marker / f"{Path(input_reads).stem}.match.fq.gz"

        with (
            pysam.FastxFile(input_reads) as fin,
            gzip.open(outfile, mode="wt") as fout,
        ):
            for read in fin:
                if read.name in aligned_reads:
                    fout.write(str(read) + "\n")

        return outfile

    def _search(
        self,
        query_db: Path,
        target_db: Path,
        outdb: str,
        min_seq_id: float = 0.9,
        min_aln_len: float = 50,
    ) -> MMSeqsEasySearchParser:
        """Search sequences against database."""
        search_result_db = self.outdir / self.marker / outdb / "db"

        search = MMSeqsEasySearchConfig(
            query_db=query_db,
            target_db=target_db,
            result_db=search_result_db,
            tmpdir=self.outdir / self.marker / "mmseqs_tmp",
            search_type=3,
            v=1,
            threads=self.threads,
            min_seq_id=min_seq_id,
            min_aln_len=min_aln_len,
        )
        search.run()
        return MMSeqsEasySearchParser(search)

    def _assemble_reads(self, reads: Path) -> Path:
        """Assemble reads and create database. Returns path to assembly DB."""
        spades_outdir = self.outdir / self.marker / "spades.out"
        spades_db = spades_outdir / "mmseqs" / "db"

        assembler = SpadesRunner(threads=self.threads, rna=self.rna)
        contigs = assembler.assemble(reads, spades_outdir)

        db_config = CreateMMSeqsDBConfig(
            fasta_file=contigs,
            sequence_db=spades_db,
            db_type=2,
            v=1,
        )
        db_config.run()

        return spades_db

    def analyse_short_reads(
        self,
        input_reads: Path,
        marker: str,
        db: Path,
        min_seq_id: float = 0.96,
        min_aln_len: int = 100,
    ):
        """Run complete analysis pipeline."""
        self.marker = marker

        # Search reads and extract matches
        reads_search = self._search(
            query_db=input_reads, target_db=db, outdb="search_reads"
        )
        self.aligned_reads = self._extract_reads(reads_search, input_reads)

        # Assemble reads and search contigs
        self.assembly_db = self._assemble_reads(self.aligned_reads)
        contigs_search = self._search(
            query_db=db,
            target_db=self.assembly_db,
            outdb="search_contigs",
            min_seq_id=min_seq_id,
            min_aln_len=min_aln_len,
        )

        # Save final results
        output_file = self.outdir / f"{marker}.contigs.search.tsv"
        shutil.copy(contigs_search, output_file)

        return output_file
