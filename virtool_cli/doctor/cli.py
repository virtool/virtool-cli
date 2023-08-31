from pathlib import Path
import click
# import structlog

from virtool_cli.doctor.checkup import run as run_checkup
from virtool_cli.doctor.fix_reference import run as run_fix_all
from virtool_cli.doctor.fix_otu import run as run_fix_otu
from virtool_cli.doctor.repair import run as run_legacy_repair
from virtool_cli.doctor.taxid import run as run_legacy_taxid

ERROR_MESSAGE = click.style("ERROR: ", fg='red')

@click.group("doc")
def doc():
    """
    Commands related to reference repairs
    """
    pass

@doc.command()
@click.option(
    "-src", "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option('--debug/--no-debug', default=False)
def checkup(src_path, debug):
    """Diagnose outstanding issues in reference. Read-only."""
    try:
        run_fix_all(Path(src_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)

@doc.command()
@click.option(
    "-src", "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option('--debug/--no-debug', default=False)
def reference(src_path, debug):
    """Attempt a repair of all data in the src directory."""
    try:
        run_fix_all(Path(src_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)

@doc.command()
@click.option(
    "-otu", "--otu_path",
    required=True,
    type=str,
    help="the path to a single OTU directory",
)
@click.option(
    "-src", "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option('--debug/--no-debug', default=False)
def otu(otu_path, src_path, debug):
    """Attempt a repair of a given OTU directory."""

    try:
        run_fix_otu(Path(otu_path), Path(src_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)

@doc.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to the input reference.json file",
)
@click.option("-f", "--force_update", is_flag=True)
@click.option('--debug/--no-debug', default=False)
def taxid(src_path, force_update, debug):
    """Fetch taxid for all OTU in given reference directory. Formerly `virtool ref taxid`"""
    try:
        click.echo("Scanning OTUs for taxon IDs...")
        run_legacy_taxid(Path(src_path), force_update, debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + "Not a valid reference directory")
        click.echo(e)

@doc.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option('--debug/--no-debug', default=False)
def legacy(src_path, debug):
    """Formerly `virtool ref repair`."""
    if not Path(src_path).exists():
        click.echo('Source directory does not exist')

    try:
        run_legacy_repair(Path(src_path), debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + "Not a valid reference directory")
        click.echo(e)


if __name__ == "__main__":
    doc()