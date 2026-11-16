import time
from datetime import timedelta
from pathlib import Path

import click
from loguru import logger

from markadoros.database.bold_tsv_processor import BoldTSVProcessor
from markadoros.database.database_creator import DatabaseCreator
from markadoros.database.header_processor import process_bold_header
from markadoros.utils.utils import set_mmseqs_path


@click.command("bold-coi-from-tsv")
@click.option(
    "--prefix",
    type=str,
    help="Prefix to prefix output database names with.",
    required=True,
)
@click.option(
    "--min-length",
    type=int,
    default=200,
    help="Minimum length of a sequence to retain.",
)
@click.option(
    "--create-index/--no-create-index",
    is_flag=True,
    default=False,
    help="Create MMSeqs2 indexes for each marker database.",
)
@click.option(
    "--exclude-file",
    type=click.Path(exists=True),
    help="New-line separated file of regular expressions to exclude header matches from the output database.",
)
@click.option(
    "--outdir",
    "-o",
    type=click.Path(),
    default=Path.cwd() / "markadoros.db",
    help="Output directory for databases. Defaults to current working directory.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=1,
    help="Number of threads to use for MMSeqs2 database creation.",
)
@click.option(
    "--cleanup/--no-cleanup",
    is_flag=True,
    default=True,
    help="Clean up temporary files after database creation.",
)
@click.argument("tsv", type=click.Path(exists=True))
def bold_tsv_database(
    prefix: str,
    min_length: int,
    tsv: str,
    create_index: bool,
    outdir: str,
    threads: int,
    cleanup: bool,
    exclude_file: str | None,
):
    """
    Build COI MMSeqs2 database from a BOLD TSV file, and record their parameters.

    TSV: Path to the TSV file to build the database from.
    """
    start_time = time.perf_counter()
    set_mmseqs_path()

    logger.add(
        str(Path(outdir).resolve() / "markadoros.database.log"),
        format="[{time:HH:mm:ss}] | markadoros | {level} - {message}",
        level="INFO",
    )

    tsv_processor = BoldTSVProcessor(
        tmpdir=Path(outdir) / "tmp",
    )
    bins_fasta, bins_taxa = tsv_processor.process_to_fasta(Path(tsv))

    database_creator = DatabaseCreator(
        outdir=Path(outdir),
        header_processor=process_bold_header,
        deduplicate=False,
        cluster=False,
        create_index=create_index,
        min_length=min_length,
        threads=threads,
        exclude_file=Path(exclude_file) if exclude_file else None,
    )
    database_creator.create_marker_database(
        fasta=bins_fasta,
        markers=["COI"],
        prefix=prefix,
        bins_taxa=bins_taxa,
    )

    elapsed = time.perf_counter() - start_time
    formatted_time = str(timedelta(seconds=int(elapsed)))

    logger.info(f"Finished! Database creation completed in {formatted_time}.")

    if cleanup:
        database_creator.cleanup()
