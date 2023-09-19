from pathlib import Path
import click

from virtool_cli.doctor.checkup import run as run_checkup
from virtool_cli.doctor.fix_reference import run as run_fix_all
from virtool_cli.doctor.fix_otu import run as run_fix_otu

ERROR_MESSAGE = click.style("ERROR: ", fg="red")


@click.group("doctor")
def doctor():
    """
    Commands related to reference repairs
    """
    pass


@doctor.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option("--debug/--no-debug", default=False)
def checkup(src_path, debug):
    """Diagnose outstanding issues in reference. Read-only."""
    try:
        run_checkup(src_path, debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


@doctor.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option("--debug/--no-debug", default=False)
def reference(src_path, debug):
    """Attempt a repair of all data in the src directory."""
    try:
        run_fix_all(src_path, debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


@doctor.command()
@click.option(
    "-otu",
    "--otu_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a single OTU directory",
)
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option("--debug/--no-debug", default=False)
def otu(otu_path, src_path, debug):
    """Attempt a repair of a given OTU directory."""

    try:
        run_fix_otu(otu_path, src_path, debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


if __name__ == "__main__":
    doctor()
