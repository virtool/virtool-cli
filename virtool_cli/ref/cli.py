import glob
import sys
from pathlib import Path

import click
from structlog import get_logger

from virtool_cli.add.cli import add
from virtool_cli.legacy.utils import iter_legacy_otus
from virtool_cli.legacy.validate import validate_legacy_repo
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.options import debug_option, path_option
from virtool_cli.ref.build import build_json
from virtool_cli.ref.repo import EventSourcedRepo
from virtool_cli.ref.resources import DataType
from virtool_cli.ref.utils import format_json
from virtool_cli.update.cli import update
from virtool_cli.utils.logging import configure_logger

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

    repo.create_otu("", None, taxonomy.name, None, [], taxid)

    logger.info("Created OTU", id=str(otu.id), name=otu.name, taxid=taxid)


@ref.group()
def legacy():
    """Validate and convert legacy references."""


@legacy.command()
@click.option(
    "--path",
    help="the name to the legacy reference repository",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def precache(path: Path):
    """Pre-cache all accessions in a legacy reference repository."""
    ncbi = NCBIClient(path / ".migration_cache", False)

    buffer = []

    for otu_ in iter_legacy_otus(path / "src"):
        for isolate in otu_["isolates"]:
            for sequence in isolate["sequences"]:
                buffer.append(sequence["accession"])

                if len(buffer) > 450:
                    ncbi.fetch_genbank_records(buffer)
                    buffer = []

    if buffer:
        ncbi.fetch_genbank_records(buffer)


@legacy.command(name="format")
@click.option(
    "--path",
    help="the name to the legacy reference repository",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def reformat(path: Path):
    """Format a legacy reference repository.

    Re-formats every JSON file in a legacy reference repository to have a consistent
    format.
    """
    src_path = path / "src"

    for json_path in glob.glob(str(src_path / "**/*.json"), recursive=True):
        format_json(Path(json_path))


@legacy.command()
@debug_option
@click.option(
    "--fix",
    is_flag=True,
    help="attempt to fix errors in-place",
)
@click.option(
    "--limit",
    default=0,
    help="exit if this many otus have errors",
)
@click.option(
    "--no-ok",
    is_flag=True,
    help="don't print anything if an otu is valid",
)
@click.option(
    "--path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a legacy reference directory",
)
def validate(debug: bool, fix: bool, limit: int, no_ok: bool, path: Path):
    """Validate a legacy reference repository."""
    configure_logger(debug)

    validate_legacy_repo(fix, limit, no_ok, path)


@ref.command()
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
@click.option(
    "-o",
    "--output-path",
    required=True,
    type=click.Path(exists=False, file_okay=True, path_type=Path),
    help="the path to write the reference.json file to",
)
@click.option(
    "-p",
    "--path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to the reference repository",
)
def build(output_path: Path, path: Path, indent: bool, version: str):
    """Build a Virtool reference.json file from a reference repository."""
    build_json(indent, output_path, path, version)
