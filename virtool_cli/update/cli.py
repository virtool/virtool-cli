from pathlib import Path
import click

from virtool_cli.update.update_ref import run as run_update_all
from virtool_cli.update.update_otu import run as run_update_single


@click.group("update")
def update():
    """
    Commands related to updates
    """
    pass


@update.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option(
    "--filter",
    default="*",
    type=str,
    help="entry point for glob filter",
)
@click.option(
    "--dry/--not-dry",
    default=False,
    help="Write to default cache directory instead of reference directory",
)
@click.option("--evaluate/--no-evaluate", default=False, help="Enable auto-evaluation")
@click.option("--debug/--no-debug", default=False, help="Enable debugging logs")
def reference(src_path, filter, evaluate, dry, debug):
    """Fetch new NCBI data for all OTU in a given reference directory."""
    try:
        run_update_all(
            src_path, filter, auto_evaluate=evaluate, dry_run=dry, debugging=debug
        )
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory", err=True)


@update.command()
@click.option(
    "-otu",
    "--otu_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a single OTU directory",
)
@click.option(
    "--dry/--not-dry",
    default=False,
    help="Write to default cache directory instead of reference directory",
)
@click.option("--evaluate/--no-evaluate", default=False, help="Enable auto-evaluation")
@click.option("--debug/--no-debug", default=False, help="Enable debugging logs")
def otu(otu_path, evaluate, dry, debug):
    """Fetch new NCBI data for a given OTU directory."""

    try:
        run_update_single(
            otu_path, auto_evaluate=evaluate, dry_run=dry, debugging=debug
        )

    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory", err=True)


if __name__ == "__main__":
    update()
