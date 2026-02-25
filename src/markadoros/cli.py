import json
import time
from datetime import timedelta
from pathlib import Path

import click

from markadoros.data import BOLD_CONFIG, UNITE_CONFIG
from markadoros.database_creator import DatabaseCreator
from markadoros.header_processor import (
    process_bold_header,
    process_unite_header,
)
from markadoros.search_pipeline import SearchPipeline
from markadoros.utils import validate_and_load_index


@click.group()
@click.version_option(
    package_name="markadoros",
    message="%(prog)s %(version)s",
)
def cli():
    pass


@cli.command()
@click.option(
    "--preset",
    "-p",
    type=click.Choice(["bold", "unite"]),
    help="Use a preset profile to generate databases.",
    default=None,
)
@click.option(
    "--db-params",
    type=click.Path(exists=True),
    help="If building a custom database, a JSON file describing each sub-database to be generated.",
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
    preset: str,
    db_params: str,
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

    if preset == "bold":
        db_dict = BOLD_CONFIG
        header_processor = process_bold_header
    elif preset == "unite":
        db_dict = UNITE_CONFIG
        header_processor = process_unite_header
    else:
        with open(db_params, "r") as f:
            db_dict = json.load(f)
        header_processor = None

    database_creator = DatabaseCreator(
        outdir=Path(outdir),
        db_dict=db_dict,
        header_processor=header_processor,
    )
    database_creator.create_marker_database(
        fasta=Path(fasta),
    )

    elapsed = time.perf_counter() - start_time
    formatted_time = str(timedelta(seconds=int(elapsed)))

    click.echo(f"Finished! Database creation completed in {formatted_time}.")

    if cleanup:
        database_creator.cleanup()


@cli.command()
@click.option(
    "--type",
    "-x",
    type=click.Choice(
        [
            "sr",
            "short",
            "illumina",
            "rnaseq",
            "pb",
            "pacbio_hifi",
            "ont",
            "oxford_nanopore",
            "contigs",
        ]
    ),
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
    default=Path.cwd(),
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
    "--include-lineage",
    is_flag=True,
    default=False,
    help="Include lineage information in output.",
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
    include_lineage: bool,
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

    # Run pipeline
    pipeline = SearchPipeline(
        outdir=Path(outdir),
        threads=threads,
        type=type,
        database_index=database_index,
        include_lineage=include_lineage,
    )
    pipeline.run(
        input=Path(input),
        n_reads=nreads,
        db_name=db,
        prefix=prefix,
    )

    elapsed = time.perf_counter() - start_time
    formatted_time = str(timedelta(seconds=int(elapsed)))

    click.echo(f"Finished! Sequences searched in {formatted_time}.")

    if cleanup:
        pipeline.cleanup()


if __name__ == "__main__":
    cli()
