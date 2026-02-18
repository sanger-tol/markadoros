import json
import time
from datetime import timedelta
from pathlib import Path

import click

from markadoros.data import BOLD_CONFIG, UNITE_CONFIG
from markadoros.database_creator import DatabaseCreator
from markadoros.search_pipeline import SearchPipeline
from markadoros.utils import validate_and_load_index


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--preset",
    type=str,
    help="Use a preset profile to generate databases.",
)
@click.option(
    "--db-params",
    type=click.Path(exists=True),
    help="If the input FASTA file contains sequences from multiple barcode genes, the path to a JSON file describing each sub-database to be generated.",
)
@click.option(
    "--outdir",
    type=click.Path(),
    default=Path.cwd() / "markadoros.db",
    help="Output directory for databases. Defaults to current working directory.",
)
@click.option(
    "--threads",
    type=int,
    default=1,
    help="Number of threads to use for MMSeqs2 database creation.",
)
@click.option(
    "--cleanup",
    is_flag=True,
    default=True,
    help="Clean up temporary files after database creation.",
)
@click.argument("fasta", type=click.Path(exists=True))
def build_database(
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
    elif preset == "unite":
        db_dict = UNITE_CONFIG
    else:
        db_dict = json.load(db_params)

    database_creator = DatabaseCreator(
        outdir=outdir,
        db_dict=db_dict,
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
    "--platform",
    type=str,
    help="Input read platform. One of illumina, illumina_rnaseq, pacbio_hifi, oxford_nanopore.",
)
@click.option(
    "--index",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the marker database index JSON file.",
)
@click.option(
    "--outdir",
    type=click.Path(exists=False),
    default=Path.cwd(),
)
@click.option(
    "--prefix",
    type=str,
    help="Output file prefix.",
)
@click.option(
    "--nreads",
    type=int,
    help="Number of reads to process. If unspecified, all reads are processed.",
)
@click.option(
    "--threads",
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
    "--cleanup",
    is_flag=True,
    default=True,
    help="Clean up temporary files after completion.",
)
@click.argument(
    "reads",
    type=click.Path(exists=True),
)
def search_reads(
    platform: str,
    index: str,
    outdir: str,
    prefix: str,
    nreads: int,
    threads: int,
    db: str,
    cleanup: bool,
    reads: str,
):
    """Search a set of short reads against a marker database.

    READS: Path to the reads file (FASTX or CRAM) to search for barcodes in.
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
        platform=platform,
        cleanup=cleanup,
        database_index=database_index,
    )

    pipeline.run(
        reads=Path(reads),
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
