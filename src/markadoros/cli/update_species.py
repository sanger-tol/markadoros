from typing import IO

import click
from loguru import logger

from markadoros.summary.results_models import ResultsReport
from markadoros.taxonomy.taxonomy import BinsTaxa, GOATSynonymGetter


@click.command("update-species")
@click.option(
    "--expected-taxon",
    type=str,
    help="The expected taxon for the search results.",
    required=True,
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
    "--output",
    "-o",
    help="Path to write the updated summary JSON file to, or stdout if not provided.",
    type=click.File("w"),
    default="-",
)
@click.argument(
    "summary_json",
    type=click.File("r"),
    help="The summary JSON file to update.",
)
def update_species(
    expected_taxon: str,
    find_goat_synonyms: bool,
    synonyms: str,
    summary_json: IO,
    output: IO,
):
    """
    Read a summary JSON file from a previous markadoros run and update the expected taxon to a new species
    """
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

    summary = ResultsReport.model_validate_json(summary_json.read())

    bins_taxa = None
    if (summary.input.database.parent / "bins_taxa.json.gz").exists():
        bins_taxa = BinsTaxa(summary.input.database.parent / "bins_taxa.json.gz")

    logger.info(
        f"Updating expexted taxon to {expected_taxon} with synonyms {input_synonyms}."
    )
    summary.update_expected_species(
        expected_taxon=expected_taxon,
        synonyms=input_synonyms,
        bin_taxa=bins_taxa,
    )

    logger.info("Writing updated summary.")
    with output as f:
        f.write(summary.model_dump_json(indent=4))
