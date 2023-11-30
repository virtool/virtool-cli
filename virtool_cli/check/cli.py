from pathlib import Path
import click

from virtool_cli.check.check_otu import run as run_otu
from virtool_cli.check.check_reference import run as run_reference

ERROR_MESSAGE = click.style("ERROR: ", fg="red")


@click.group("check")
def check():
    """
    Commands related to validation checks
    """
    pass


@check.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option("--debug/--no-debug", default=False)
def reference(src_path, debug):
    """Checks the validity of a reference source"""
    try:
        run_reference(src_path, debug)

    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            ERROR_MESSAGE + f"{str(src_path)} is not a valid reference directory",
            err=True,
        )


@check.command()
@click.option(
    "-otu",
    "--otu_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a single OTU directory",
)
@click.option("--debug/--no-debug", default=False)
def otu(otu_path, debug):
    """Checks the validity of a reference source"""
    try:
        run_otu(otu_path, debug)

    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            ERROR_MESSAGE + f"{str(otu_path)} is not a valid OTU directory",
            err=True,
        )


if __name__ == "__main__":
    check()
