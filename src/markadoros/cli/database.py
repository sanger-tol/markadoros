import time
from datetime import timedelta
from pathlib import Path

import click
from loguru import logger

from markadoros.database.database_creator import DatabaseCreator
from markadoros.database.header_processor import (
    process_bold_header,
    process_generic_header,
    process_silva_lsu_header,
    process_silva_ssu_header,
    process_unite_header,
)
from markadoros.utils.utils import set_mmseqs_path


@click.command("database")
@click.option(
    "--header-type",
    "-x",
    type=click.Choice(["bold", "unite", "silva_lsu", "silva_ssu", "generic"]),
    help="Use a preset header processor to generate databases.",
    default="generic",
)
@click.option(
    "--markers",
    type=str,
    help="Comma-separated list of marker gene names. Required if --header-type is `generic`, otherwise only used if `--header_type` is `bold`.",
    required=False,
)
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
    "--deduplicate/--no-deduplicate",
    is_flag=True,
    default=False,
    help="Whether to deduplicate sequences in the output database.",
)
@click.option(
    "--cluster/--no-cluster",
    is_flag=True,
    default=False,
    help="Cluster the database sequences with mmseqs linclust.",
)
@click.option(
    "--cluster-min-seq-id",
    type=float,
    default=0.99,
    help="Percent identity at which to cluster sequences",
)
@click.option(
    "--cluster-coverage",
    type=float,
    default=0.8,
    help="Coverage overlap at which to cluster sequences",
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
@click.argument("fasta", type=click.Path(exists=True))
def database(
    header_type: str,
    markers: str | None,
    prefix: str,
    min_length: int,
    deduplicate: bool,
    cluster: bool,
    cluster_coverage: float,
    cluster_min_seq_id: float,
    fasta: str,
    create_index: bool,
    outdir: str,
    threads: int,
    cleanup: bool,
    exclude_file: str | None,
):
    """
    Build MMSeqs2 databases from a FASTA file, and record their parameters.

    FASTA: Path to the FASTA file to build the database from.
    """
    start_time = time.perf_counter()
    set_mmseqs_path()

    logger.add(
        str(Path(outdir).resolve() / "markadoros.database.log"),
        format="[{time:HH:mm:ss}] | markadoros | {level} - {message}",
        level="INFO",
    )

    if header_type == "bold":
        header_processor = process_bold_header
        if markers is None:
            markers_list = ["COI", "CYTB", "rbcL", "matK", "18S", "28S"]
        else:
            markers_list = markers.split(",")
    elif header_type == "unite":
        header_processor = process_unite_header
        markers_list = ["ITS"]
    elif header_type == "silva_lsu":
        header_processor = process_silva_lsu_header
        markers_list = ["LSU"]
    elif header_type == "silva_ssu":
        header_processor = process_silva_ssu_header
        markers_list = ["SSU"]
    else:
        header_processor = process_generic_header

    if markers is None:
        logger.info(
            f"--header-type is {header_type}! Setting marker to {', '.join(markers_list)}!"
        )

    database_creator = DatabaseCreator(
        outdir=Path(outdir),
        header_processor=header_processor,
        deduplicate=deduplicate,
        cluster=cluster,
        cluster_coverage=cluster_coverage,
        cluster_min_seq_id=cluster_min_seq_id,
        create_index=create_index,
        min_length=min_length,
        threads=threads,
        exclude_file=Path(exclude_file) if exclude_file else None,
    )
    database_creator.create_marker_database(
        fasta=Path(fasta),
        markers=markers_list,
        prefix=prefix,
    )

    elapsed = time.perf_counter() - start_time
    formatted_time = str(timedelta(seconds=int(elapsed)))

    logger.info(f"Finished! Database creation completed in {formatted_time}.")

    if cleanup:
        database_creator.cleanup()
