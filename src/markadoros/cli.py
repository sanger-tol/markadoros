import json
import sys
import time
from datetime import timedelta
from pathlib import Path

import click
from loguru import logger

from markadoros.database_creator import DatabaseCreator
from markadoros.header_processor import (
    process_bold_header,
    process_generic_header,
    process_unite_header,
)
from markadoros.input_types import get_valid_input_types, normalize_input_type
from markadoros.search_pipeline import SearchPipeline
from markadoros.utils import set_mmseqs_path, validate_and_load_index


@click.group()
@click.version_option(
    package_name="markadoros",
    message="%(prog)s %(version)s",
)
def cli():
    logger.remove()
    logger.add(
        sys.stderr,
        format="[{time:HH:mm:ss}] | markadoros | {level} - {message}",
        level="INFO",
    )
    set_mmseqs_path()


@cli.command()
@click.option(
    "--header-type",
    "-x",
    type=click.Choice(["bold", "unite"]),
    help="Use a preset header processor to generate databases.",
    default=None,
)
@click.option(
    "--marker",
    type=str,
    help="The name of a marker gene to find in the input FASTA. Can be specified multiple times for each marker present in the data.",
    multiple=True,
    required=True,
)
@click.option(
    "--prefix",
    type=str,
    help="Prefix to prefix output database names with.",
    required=True,
)
@click.option(
    "--min_length",
    type=int,
    default=200,
    help="Minimum length of a sequence to retain.",
)
@click.option(
    "--deduplicate/--no-deduplicate",
    is_flag=True,
    default=True,
    help="Whether to deduplicate sequences in the output database.",
)
@click.option(
    "--cluster/--no-cluster",
    is_flag=True,
    default=True,
    help="Cluster the database sequences with mmseqs linclust.",
)
@click.option(
    "--cluster_min_seq_id",
    type=float,
    default=0.99,
    help="Percent identity at which to cluster sequences",
)
@click.option(
    "--cluster_coverage",
    type=float,
    default=0.8,
    help="Coverage overlap at which to cluster sequences",
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
@click.argument("fasta", type=click.Path(exists=True))
def database(
    header_type: str,
    marker: list[str],
    prefix: str,
    min_length: int,
    deduplicate: bool,
    cluster: bool,
    cluster_coverage: float,
    cluster_min_seq_id: float,
    fasta: str,
    outdir: str,
    threads: int,
    cleanup: bool,
):
    """
    Build MMSeqs2 databases from a FASTA file, and record their parameters.

    FASTA: Path to the FASTA file to build the database from.
    """
    start_time = time.perf_counter()

    logger.add(
        str(Path(outdir).resolve() / "markadoros.database.log"),
        format="[{time:HH:mm:ss}] | markadoros | {level} - {message}",
        level="INFO",
    )

    if header_type == "bold":
        header_processor = process_bold_header
    elif header_type == "unite":
        header_processor = process_unite_header
        logger.warning(
            "--header-type is unite! All --marker specifications are ignored and 'ITS' is forced!"
        )
        marker = ["ITS"]
    else:
        header_processor = process_generic_header

    database_creator = DatabaseCreator(
        outdir=Path(outdir),
        header_processor=header_processor,
        deduplicate=deduplicate,
        cluster=cluster,
        cluster_coverage=cluster_coverage,
        cluster_min_seq_id=cluster_min_seq_id,
        min_length=min_length,
        threads=threads,
    )
    database_creator.create_marker_database(
        fasta=Path(fasta),
        markers=list(marker),
        prefix=prefix,
    )

    elapsed = time.perf_counter() - start_time
    formatted_time = str(timedelta(seconds=int(elapsed)))

    logger.info(f"Finished! Database creation completed in {formatted_time}.")

    if cleanup:
        database_creator.cleanup()


@cli.command()
@click.option(
    "--type",
    "-x",
    type=click.Choice(get_valid_input_types()),
    default="sr",
    required=True,
    help="Input data type, either a sequencing platform (to choose the right assembler), or pre-assembled contigs.",
)
@click.option(
    "--index",
    "-i",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the marker database index JSON file.",
)
@click.option(
    "--outdir",
    "-o",
    type=click.Path(exists=False),
    default=str(Path.cwd().resolve()),
)
@click.option(
    "--prefix",
    "-p",
    type=str,
    help="Output file prefix.",
)
@click.option(
    "--nreads",
    "-n",
    type=int,
    help="Number of reads to process. If unspecified, all reads are processed.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=1,
    help="Number of threads to use for searching and assembly.",
)
@click.option(
    "--db",
    type=str,
    help="Optionally, the name of a single database within the index file to search with.",
)
@click.option(
    "--cleanup/--no-cleanup",
    is_flag=True,
    default=True,
    help="Clean up temporary files after completion.",
)
@click.option(
    "--min_seq_id",
    "-m",
    type=float,
    default=0.96,
    help="Minimum sequence ID required to report a hit.",
)
@click.option(
    "--min_aln_len",
    "-l",
    type=int,
    default=450,
    help="Minimum alignment length required to report a hit.",
)
@click.option(
    "--expected_taxon",
    type=str,
    help="The expected taxon binomial name.",
)
@click.argument(
    "input",
    type=click.Path(exists=True),
)
def search(
    type: str,
    index: str,
    outdir: str,
    prefix: str,
    nreads: int,
    expected_taxon: str,
    min_seq_id: float,
    min_aln_len: int,
    threads: int,
    db: str,
    cleanup: bool,
    input: str,
):
    """Search reads or contigs against a marker database.

    INPUT: Path to reads file (FASTX or CRAM) or contigs file (FASTA).
    """
    start_time = time.perf_counter()

    # Load and validate database
    try:
        database_index = validate_and_load_index(Path(index))
    except (FileNotFoundError, json.JSONDecodeError) as e:
        raise click.ClickException(f"Could not load database from {index}: {e}")

    # Normalize input type
    input_type = normalize_input_type(type)

    # Run pipeline
    pipeline = SearchPipeline(
        outdir=Path(outdir),
        threads=threads,
        input_type=input_type,
        database_index=database_index,
        expected_taxon=expected_taxon,
        min_seq_id=min_seq_id,
        min_aln_len=min_aln_len,
    )
    pipeline.run(
        input=Path(input),
        n_reads=nreads,
        db_name=db,
        prefix=prefix,
    )

    elapsed = time.perf_counter() - start_time
    formatted_time = str(timedelta(seconds=int(elapsed)))

    logger.info(f"Finished! Sequences searched in {formatted_time}.")

    if cleanup:
        pipeline.cleanup()


if __name__ == "__main__":
    cli()
