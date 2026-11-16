import time
from datetime import timedelta
from pathlib import Path

import click
from loguru import logger

from markadoros.database.database_index import DatabaseIndex
from markadoros.search.search_pipeline import SearchPipeline
from markadoros.summary.results_models import ExpectedTaxonInfo
from markadoros.taxonomy.taxonomy import GOATSynonymGetter
from markadoros.utils.input_types import get_valid_input_types, normalize_input_type
from markadoros.utils.utils import (
    get_simple_name,
    set_mmseqs_path,
)


@click.command()
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
    required=True,
)
@click.option(
    "--cleanup/--no-cleanup",
    is_flag=True,
    default=True,
    help="Clean up temporary files after completion.",
)
@click.option(
    "--db-to-tmpdir/--no-db-to-tmpdir",
    is_flag=True,
    default=True,
    help="Copy the database to the tempdir. Can be useful if running multiple searches simultaneously on the same database, which can cause IO issues.",
)
@click.option(
    "--min-seq-id",
    "-m",
    type=float,
    default=0.96,
    help="Minimum sequence ID required to report a hit.",
)
@click.option(
    "--min-aln-len",
    "-l",
    type=int,
    default=250,
    help="Minimum alignment length required to report a hit.",
)
@click.option(
    "--expected-taxon",
    "-e",
    type=str,
    help="The expected taxon binomial name.",
)
@click.option(
    "--find-goat-synonyms",
    "-g",
    is_flag=True,
    default=False,
    help="Find synonyms for the expected taxon on GoaT.",
)
@click.option(
    "--synonyms",
    "-s",
    help="A comma-separated list of synonyms to check in addition to the expected taxon. If provided, GoaT lookup is disabled.",
)
@click.option(
    "--save-contigs",
    is_flag=True,
    default=False,
    help="Save the assembled contigs to file in the output directory.",
)
@click.option(
    "--fallback-reads",
    type=int,
    default=30,
    help="If assembly fails, fall back and search the first N matched reads for a match.",
)
@click.argument(
    "input",
    type=click.Path(exists=True),
)
def search(
    type: str,
    index: str,
    outdir: str,
    prefix: str | None,
    nreads: int,
    expected_taxon: str,
    find_goat_synonyms: bool,
    synonyms: str,
    min_seq_id: float,
    min_aln_len: int,
    threads: int,
    db: str,
    db_to_tmpdir: bool,
    cleanup: bool,
    save_contigs: bool,
    input: str,
    fallback_reads: int,
):
    """Search reads or contigs against a marker database.

    INPUT: Path to reads file (FASTX or CRAM) or contigs file (FASTA).
    """
    start_time = time.perf_counter()
    set_mmseqs_path()

    # Validate fallback reads
    if fallback_reads < 0:
        raise click.ClickException("--fallback-reads must be greater than 0.")

    # Load and validate database
    try:
        database_index = DatabaseIndex(Path(index)).load()
    except (FileNotFoundError, ValueError) as e:
        raise click.ClickException(f"Could not load database from {index}: {e}")

    # Handle synonym finding
    if synonyms:
        if find_goat_synonyms:
            logger.warning(
                "Both --find-goat-synonyms and --synonyms provided; using --synonyms only."
            )
        input_synonyms = synonyms.split(",")
    elif find_goat_synonyms:
        goat_lookup = GOATSynonymGetter()
        input_synonyms = goat_lookup.get_synonyms(expected_taxon)
    else:
        input_synonyms = []

    # Normalize input type
    input_type = normalize_input_type(type)

    # Initialise prefix
    if prefix is None:
        prefix = get_simple_name(Path(input))

    tmpdir = Path(outdir) / f"{prefix}.tmp"

    expected = (
        ExpectedTaxonInfo(taxon=expected_taxon, synonyms=input_synonyms)
        if expected_taxon
        else None
    )

    # Run pipeline
    pipeline = SearchPipeline(
        outdir=Path(outdir),
        tmpdir=tmpdir,
        threads=threads,
        input_type=input_type,
        database_index=database_index,
        db_to_tmpdir=db_to_tmpdir,
        expected_taxon=expected,
        min_seq_id=min_seq_id,
        min_aln_len=min_aln_len,
        save_contigs=save_contigs,
        fallback_reads=fallback_reads,
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
