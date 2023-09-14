from pathlib import Path
import click

from virtool_cli.update.update_ref import run as run_update_all
from virtool_cli.update.update import run as run_update_single
from virtool_cli.update.isolate import run as run_isolate

@click.group("update")
def update():
    """
    Commands related to updates
    """
    pass

@update.command()
@click.option(
    "-src", "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=str,
    default='.cache/catalog',
    help="the path to a catalog directory",
)
@click.option('--filter/--no-filter', default=False)
@click.option('--debug/--no-debug', default=False)
def reference(src_path, catalog_path, filter, debug):
    """Fetch new sequences and isolates for all OTU in a given reference directory."""
    if not Path(catalog_path).exists():
        click.echo("Not a valid catalog directory")
        return

    try:
        run_update_all(Path(src_path), Path(catalog_path), auto_evaluate=filter, debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)

@update.command()
@click.option(
    "-otu", "--otu_path",
    required=True,
    type=str,
    help="the path to a single OTU directory",
)
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=str,
    default='.cache/catalog',
    help="the path to a catalog directory",
)
@click.option('--filter/--no-filter', default=False)
@click.option('--debug/--no-debug', default=False)
def otu(otu_path, catalog_path, filter, debug):
    """Fetch new sequences and isolates for a given OTU directory."""
    if not Path(catalog_path).exists():
        click.echo("Not a valid catalog directory")
        return

    try:
        run_update_single(
            Path(otu_path), Path(catalog_path), 
            auto_evaluate=filter, debugging=debug
        )

    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)

@update.command()
@click.option(
    "-src", "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option('--debug/--no-debug', default=False)
def legacy(src_path, debug):
    """Formerly `virtool ref isolate`. Use for quick retrival."""
    try:
        run_isolate(Path(src_path))
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


if __name__ == "__main__":
    update()