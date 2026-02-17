import json
import shutil
from pathlib import Path

import click

from markadoros.data import BOLD_CONFIG, UNITE_CONFIG
from markadoros.databases import DatabaseCreator, validate_db
from markadoros.read_preprocessor import ReadPreprocessor
from markadoros.short_read_analyser import ShortReadAnalyser
from markadoros.spades import SpadesRunner


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
    "--db_params",
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
    preset,
    db_params,
    fasta,
    outdir,
    threads,
    cleanup,
):
    """
    Build MMSeqs2 databases from a FASTA file, and record their parameters.
    """
    if preset == "bold":
        db_dict = BOLD_CONFIG
    elif preset == "unite":
        db_dict = UNITE_CONFIG
    else:
        db_dict = json.load(db_params)

    fasta = Path(fasta)

    bold_processor = DatabaseCreator(
        outdir=outdir,
        db_dict=db_dict,
    )
    bold_processor.create_marker_database(
        fasta=fasta,
    )

    if cleanup:
        shutil.rmtree(outdir / "tmp")


@cli.command()
@click.option("--db", type=click.Path(exists=True), required=True)
@click.option("--outdir", type=click.Path(exists=True), default=Path.cwd())
@click.argument("long_reads", type=click.Path(exists=True))
def search_long_reads():
    pass


@cli.command()
@click.option(
    "--rna",
    is_flag=True,
    help="Input sequences for searching are RNA. SPAdes will be run in RNA assembly mode.",
)
@click.option(
    "--index",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the database JSON file.",
)
@click.option(
    "--outdir", type=click.Path(exists=False), default=Path.cwd() / "markadoros.out"
)
@click.option(
    "--prefix",
    type=str,
)
@click.option("--nreads", type=int, default=10000000)
@click.option("--threads", type=int, default=1)
@click.option("--db", type=str)
@click.option(
    "--cleanup",
    is_flag=True,
    help="Remove temporary files after processing.",
    default=True,
)
@click.argument("short_reads", type=click.Path(exists=True))
def search_short_reads(
    rna, index, outdir, prefix, nreads, threads, db, cleanup, short_reads
):
    """
    Search a set of short reads against a marker database, extract the results,
    assemble then with SPAdes and then search the assembled contigs.
    """
    ## Load the database JSON
    try:
        index = validate_db(Path(index))
    except (FileNotFoundError, json.JSONDecodeError):
        raise ValueError(f"Could not load database from {index}!")

    outdir = Path(outdir)
    short_reads = Path(short_reads)

    ## Subsample the input reads
    preprocessed_reads = ReadPreprocessor(
        input=short_reads,
        outdir=outdir / "tmp",
        n_reads=nreads,
    )
    subsampled_reads_db = preprocessed_reads.subsample_reads()

    short_read_analysis = ShortReadAnalyser(
        outdir=outdir,
        threads=threads,
        prefix=prefix if prefix is not None else short_reads.stem,
        assembler=SpadesRunner(threads=threads, rna=rna),
        rmtmp=cleanup,
    )

    ## Run across all databases unless specific one specified
    databases_list = index.keys()
    if db is not None:
        if db not in databases_list:
            raise ValueError(f"Database {db} not found in index!")
        databases_list = [db]

    results = {}
    for db, params in index.items():
        if db in databases_list:
            result = short_read_analysis.analyse_short_reads(
                input_reads=subsampled_reads_db,
                marker=params.get("marker"),
                db=Path(params.get("db")),
                min_seq_id=params.get("min_seq_id"),
                min_aln_len=params.get("min_aln_len"),
            )
            if result is not None:
                results[db] = result

    ## Print some result summary
    for db, result in results.items():
        click.echo()
        click.echo(f"Top result for {db}:")
        out = result.head(1).reset_index()
        for _, row in out.iterrows():
            click.echo(f"target: {row['target']}")
            click.echo(f"query: {row['query']}")
            click.echo(f"fident: {row['fident']}")
            click.echo(f"alnlen: {row['alnlen']}")
            click.echo(f"coverage: {row['coverage']}x")
            click.echo()

    if cleanup:
        shutil.rmtree(outdir / "tmp")


if __name__ == "__main__":
    cli()
