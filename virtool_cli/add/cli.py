import asyncio
from pathlib import Path

import click

from virtool_cli.add.accessions import add_accessions
from virtool_cli.add.otu import add_otu
from virtool_cli.options import debug_option, path_option
from virtool_cli.ref.repo import EventSourcedRepo as Repo
from virtool_cli.utils.logging import configure_logger


@click.group("add")
def add():
    """Commands related to new data"""


@add.command()
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
)
@click.option("--taxid", type=int, required=True)
@path_option
@debug_option
def accessions(debug, path, taxid, accessions_: list[str]):
    """Fetch and write the data for the given NCBI accessions to an OTU.

    virtool add accessions --taxid 2697049 MN996528.1 --path [repo_path]

    """
    configure_logger(debug)

    repo = Repo(path)

    otu_index = repo.index_otus()
    if taxid not in otu_index:
        add_otu(repo, taxid)

    add_accessions(repo, taxid, accessions_)


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
