from pathlib import Path

import click

from virtool_cli.add.cli import add
from virtool_cli.check.cli import check
from virtool_cli.options import path_option
from virtool_cli.ref.build import build_json
from virtool_cli.ref.divide import run as run_divide
from virtool_cli.ref.init import init_reference
from virtool_cli.ref.migrate import run as run_migrate
from virtool_cli.update.cli import update
from virtool_cli.utils.logging import configure_logger, error_message_style


@click.group("ref")
def ref():
    """Manage references"""


ref.add_command(update)
ref.add_command(add)
ref.add_command(check)


@ref.command()
@click.option(
    "--path",
    default=".",
    help="the path to initialize the repository at",
    type=click.Path(path_type=Path),
)
@click.option("--debug", default=False, is_flag=True)
def init(debug: bool, path: Path):
    """Instantiate directory structure for an empty reference source"""
    configure_logger(debug)
    init_reference(path)


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
        click.echo(f"{error_message_style}{file_path} is not a JSON file")

    try:
        run_divide(file_path, output_path, debug)
    except (TypeError, FileNotFoundError):
        click.echo(
            f"{error_message_style}{file_path} is not a proper JSON file",
            err=True,
        )


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
