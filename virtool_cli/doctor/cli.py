from pathlib import Path
import click
# import structlog

from virtool_cli.doctor.fix_reference import run as run_fix_all
from virtool_cli.doctor.fix_otu import run as run_fix_otu

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


if __name__ == "__main__":
    doc()