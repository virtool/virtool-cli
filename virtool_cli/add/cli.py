import asyncio
from pathlib import Path

import click

from virtool_cli.add.accessions import add_accession, add_accessions
from virtool_cli.add.otus import add_otu
from virtool_cli.options import debug_option, path_option
from virtool_cli.utils.logging import configure_logger


@click.group("add")
def add():
    """Commands related to new data"""


@add.command()
@click.argument(
    "accession_",
    metavar="ACCESSION",
    type=str,
)
@path_option
@debug_option
def accession(accession_: str, debug: bool, path: Path):
    """Retrieves data for the accession."""
    configure_logger(debug)


    try:
        asyncio.run(add_accession(accession_, path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory", err=True)


@add.command()
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
)
@click.option(
    "--otu-path",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@debug_option
def accessions(accessions_: list[str], otu_path, debug):
    """Fetch and write the data for the given NCBI accessions to an OTU.

    virtool add accessions --otu-path ./tmv MN908947.3 MN996528.1

    """
    configure_logger(debug)

    try:
        asyncio.run(add_accessions(list(accessions_), otu_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("could not find otu directory", err=True)


@add.command()
@click.argument(
    "taxid",
    type=int,
)
@debug_option
@path_option
def otu(debug: bool, path: Path, taxid: int):
    """Add a new OTU using a NCBI taxonomy ID."""
    configure_logger(debug)

    try:
        asyncio.run(add_otu(taxid, path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid OTU directory", err=True)


def split_clean_csv_string(input_string: str, delimiter: str = ",") -> list[str]:
    """Splits comma-separated values into a list of strings.

    :param input_string: A raw comma-delineated string
    :param delimiter: A delimiter string for use with split(). Defaults to a comma (",")
    :return: A list of accession strings
    """
    raw_list = input_string.split(delimiter)
    stripped_list = []
    for item in raw_list:
        stripped_list.append(item.strip())

    return stripped_list
