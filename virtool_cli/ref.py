from pathlib import Path
import click

from virtool_cli.update.cli import update
from virtool_cli.add.cli import add
from virtool_cli.check.cli import check
from virtool_cli.init import run as run_init
from virtool_cli.build import run as run_build
from virtool_cli.divide import run as run_divide
from virtool_cli.migrate import run as run_migrate

ERROR_MESSAGE = click.style("ERROR: ", fg="red")


@click.group("ref")
def ref():
    """
    Commands related to reference files.
    """
    pass


ref.add_command(update)
ref.add_command(add)
ref.add_command(check)


@ref.command()
@click.option(
    "-repo",
    "--repo_path",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="the path to a containing directory for the repository",
)
@click.option("--debug/--no-debug", default=False)
def init(repo_path, debug):
    """Instantiate directory structure for an empty reference source"""
    try:
        run_init(repo_path, debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            ERROR_MESSAGE
            + "Ran into problems with the given reference repository directory",
            err=True,
        )


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a database reference directory",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default="reference.json",
    help="the output path for a reference.json file",
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
    default=None,
    type=str,
    help="the version string to include in the reference.json file",
)
@click.option("--debug/--no-debug", default=False)
def build(src_path, output, indent, version, debug):
    """Build a Virtool reference JSON file from a source directory."""
    try:
        run_build(src_path, output, indent, version, debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            ERROR_MESSAGE + "Source directory contains critical errors", err=True
        )


@ref.command()
@click.option(
    "-f",
    "--file_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="the path to a input reference.json file",
)
@click.option(
    "-o",
    "--output_path",
    default="src",
    type=click.Path(file_okay=False, path_type=Path),
    help="the output path for a divided reference directory tree",
)
@click.option("--debug/--no-debug", default=False)
def divide(file_path, output_path, debug):
    """Divide a reference.json file from Virtool into a reference directory tree."""
    if file_path.suffix != ".json":
        click.echo(ERROR_MESSAGE + f"{file_path} is not a JSON file")

    try:
        run_divide(file_path, output_path, debug)
    except (TypeError, FileNotFoundError) as e:
        click.echo(ERROR_MESSAGE + f"{file_path} is not a proper JSON file", err=True)


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

    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(
            ERROR_MESSAGE + f"{src_path} is not a valid reference directory", err=True
        )


if __name__ == "__main__":
    ref()
