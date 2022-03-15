from pathlib import Path

import click

import virtool_cli.build
import virtool_cli.divide
import virtool_cli.isolate
import virtool_cli.repair
import virtool_cli.taxid


@click.group("ref")
def ref():
    """Commands related to references"""
    pass


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a database reference directory",
)
@click.option(
    "-o",
    "--output",
    default="reference.json",
    help="the output path for a reference.json file",
)
@click.option("-i", "--indent", is_flag=True)
@click.option(
    "-V",
    "--version",
    default=None,
    type=str,
    help="the version string to include in the reference.json file",
)
def build(src_path, output, indent, version):
    """Build a Virtool reference JSON file from a reference directory."""
    try:
        virtool_cli.build.run(Path(src_path), Path(output), indent, version)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a input reference.json file",
)
@click.option(
    "-o",
    "--output",
    default="src",
    type=str,
    help="the output path for a divided reference directory tree",
)
def divide(src_path, output):
    """Divide a reference.json file from Virtool into a reference directory tree."""
    try:
        if not src_path.endswith(".json"):
            raise TypeError
        virtool_cli.divide.run(Path(src_path), Path(output))
    except (TypeError, FileNotFoundError):
        click.echo(
            "Specified reference file either does not exist or is not a proper JSON file"
        )


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to the input reference.json file",
)
@click.option("-f", "--force_update", is_flag=True)
def taxid(src_path, force_update):
    """Fetch taxid for all OTU in given reference directory."""
    try:
        virtool_cli.taxid.run(Path(src_path), force_update)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
def isolate(src_path):
    """Fetch new isolates for all OTU in a given reference directory."""
    try:
        virtool_cli.isolate.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
def repair(src_path):
    """Fix every OTU in a given reference directory."""
    try:
        virtool_cli.repair.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


if __name__ == "__main__":
    ref()
