import sys
from pathlib import Path

import click
from structlog import get_logger

from virtool_cli.add.cli import add
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.options import debug_option, path_option
from virtool_cli.ref.build import build_json
from virtool_cli.ref.migrate import run as run_migrate
from virtool_cli.ref.repo import EventSourcedRepo
from virtool_cli.ref.utils import DataType
from virtool_cli.update.cli import update
from virtool_cli.utils.logging import configure_logger, error_message_style

logger = get_logger()


@click.group("ref")
def ref():
    """Manage references"""


ref.add_command(update)
ref.add_command(add)


@ref.command()
@click.option(
    "--data-type",
    help="the type of data the reference contains (eg. genome)",
    required=True,
    type=click.Choice(DataType),
)
@click.option(
    "--name",
    help="the type of data the reference contains (eg. genome)",
    required=True,
    type=str,
)
@click.option(
    "--organism",
    default="",
    help="the organism the reference is for (eg. virus)",
    type=str,
)
@click.option(
    "--path",
    default=".",
    help="the path to initialize the repository at",
    type=click.Path(path_type=Path),
)
def init(data_type: DataType, name: str, organism: str, path: Path):
    """Create a new event-sourced repo."""
    EventSourcedRepo.new(data_type, name, path, organism)


@ref.group()
def otu():
    """Manage OTUs"""


@otu.command()
@click.argument("TAXID", type=int)
@debug_option
@path_option
def create(debug: bool, path: Path, taxid: int):
    configure_logger(debug)

    repo = EventSourcedRepo(path)
    ncbi = NCBIClient.from_repo(repo.path, False)

    taxonomy = ncbi.fetch_taxonomy_record(taxid)

    if taxonomy is None:
        logger.fatal(f"Taxonomy ID {taxid} not found")
        sys.exit(1)

    otu = repo.create_otu(
        acronym="",
        legacy_id=None,
        name=taxonomy.name,
        molecule=None,
        schema=[],
        taxid=taxid,
    )

    logger.info("Created OTU", id=str(otu.id), name=otu.name, taxid=taxid)


@ref.command()
@click.option(
    "-o",
    "--output-path",
    type=click.Path(exists=False, path_type=Path),
    default="reference.json",
    help="the output path for the reference.json file",
)
@click.option(
    "-i",
    "--indent",
    is_flag=True,
    help="auto-indent the output JSON file",
)
@click.option(
    "-V",
    "--version",
    default="",
    type=str,
    help="a version string to include in the reference.json file",
)
@path_option
def build(output_path: Path, path: Path, indent: bool, version: str):
    """Build a Virtool reference.json file from a reference repository."""
    build_json(indent, output_path, path, version)


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option("--debug/--no-debug", default=False)
def migrate(src_path, debug):
    """Convert a reference directory from v1.x to v2.x"""
    try:
        run_migrate(Path(src_path), debug)

    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            f"{error_message_style}{src_path} is not a valid reference directory",
            err=True,
        )
