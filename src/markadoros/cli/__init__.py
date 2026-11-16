import sys

import click
from loguru import logger

from markadoros.cli import bold_tsv, database, search, update_species


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
    click.echo(
        """

    ▗▖  ▗▖▗▞▀▜▌ ▄▄▄ █  ▄ ▗▞▀▜▌▐▌ ▄▄▄   ▄▄▄ ▄▄▄   ▄▄▄
    ▐▛▚▞▜▌▝▚▄▟▌█    █▄▀  ▝▚▄▟▌▐▌█   █ █   █   █ ▀▄▄
    ▐▌  ▐▌     █    █ ▀▄   ▗▞▀▜▌▀▄▄▄▀ █   ▀▄▄▄▀ ▄▄▄▀
    ▐▌  ▐▌          █  █   ▝▚▄▟▌

    """,
        err=True,
    )


cli.add_command(database.database)
cli.add_command(bold_tsv.bold_tsv_database)
cli.add_command(search.search)
cli.add_command(update_species.update_species)
