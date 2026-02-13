import json
from pathlib import Path

import click

from markadoros.databases import DatabaseCreator, validate_db
from markadoros.read_preprocessor import ReadPreprocessor
from markadoros.short_read_analyser import ShortReadAnalyser


@click.group()
def cli():
    pass


@cli.command()
@click.option("--bold", is_flag=True, default=False)
@click.option("--marker", type=str, required=False)
@click.option("--outdir", type=click.Path(exists=True), default=Path.cwd())
@click.option("--threads", type=int, default=1)
@click.option("--min_len", type=int, default=300)
@click.argument("fasta", type=click.Path(exists=True))
def build_database(bold, marker, fasta, mmseqs_path, outdir, threads, min_len):
    """Build MMSeqs databases from a FASTA release."""

    if not marker and not bold:
        raise ValueError("Either --bold or --marker must be specified")

    if marker and bold:
        click.echo(
            f"Input FASTA is a BOLD FASTA (--bold): --marker {marker} will be ignored."
        )

    bold_processor = DatabaseCreator(
        outdir=outdir,
        min_len=min_len,
    )
    bold_processor.create_marker_database(
        fasta=fasta, markers=marker, is_bold_fasta=bold
    )


@cli.command()
@click.option("--db", type=click.Path(exists=True), required=True)
@click.option("--outdir", type=click.Path(exists=True), default=Path.cwd())
@click.argument("long_reads", type=click.Path(exists=True))
def search_long_reads(ctx):
    pass


@cli.command()
@click.option("--rna", is_flag=True)
@click.option("--db", type=click.Path(exists=True), required=True)
@click.option(
    "--outdir", type=click.Path(exists=True), default=Path.cwd() / "markadoros"
)
@click.option("--n_reads", type=int, default=10000000)
@click.option("--threads", type=int, default=1)
@click.option("--markers", type=str, default=1)
@click.option("--min_seq_id", type=float, default=0.96)
@click.option("--min_aln_len", type=int, default=250)
@click.argument("short_reads", type=click.Path(exists=True))
def search_short_reads(
    short_reads, rna, db, outdir, nreads, threads, markers, min_seq_id, min_aln_len
):
    """
    Search a set of short reads against a marker database, extract the results,
    assemble then with SPAdes and then search the assembled contigs.
    """
    ## Load the database JSON
    try:
        db = json.load(validate_db(db))
    except (FileNotFoundError, json.JSONDecodeError):
        raise ValueError(f"Could not load database from {db}!")

    ## Subsample the input reads
    preprocessed_reads = ReadPreprocessor(
        input=short_reads,
        outdir=outdir / "tmp",
        n_reads=nreads,
    )
    subsampled_reads_db = preprocessed_reads.subsample_reads()

    short_read_analysis = ShortReadAnalyser(
        outdir=outdir,
        rrna=rna,
        threads=threads,
    )

    ## Run for all markers unless markers specified
    markers_list = db.keys()
    if markers:
        markers_list = markers.split(",")

    for marker in markers_list:
        db = Path(db[marker])
        short_read_analysis.analyse_short_reads(
            input_reads=subsampled_reads_db,
            marker=marker,
            db=db,
            min_seq_id=min_seq_id,
            min_aln_len=min_aln_len,
        )


if __name__ == "__main__":
    cli()
