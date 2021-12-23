from pathlib import Path

import click

import virtool_cli.build
import virtool_cli.divide
import virtool_cli.isolate
import virtool_cli.merge
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


@ref.command()
@click.option(
    "-s",
    "--source-src-path",
    required=True,
    type=str,
    help="Path to src for reference to copy isolates from",
)
@click.option(
    "-t",
    "--target-src-path",
    required=True,
    type=str,
    help="Path to src for reference to write to",
)
@click.option(
    "-r", "--resume", is_flag=True, help="Resume in-progress merge using cache"
)
@click.option(
    "-i",
    "--in-place",
    is_flag=True,
    help="Place new isolates directly into target reference",
)
@click.option(
    "-o", "--output", default="merged_ref", help="Name for merged reference directory"
)
@click.option(
    "--output-unknown-isolates",
    is_flag=True,
    help="Write CSV file for isolates with unknown source name or source type",
)
@click.option(
    "--unknown-output",
    default="unknown_isolates.txt",
    type=str,
    help="Name for CSV file containing all "
    "isolates with unknown source name or "
    "source type",
)
def merge(
    source_src_path,
    target_src_path,
    resume,
    in_place,
    output,
    output_unknown_isolates,
    unknown_output,
):
    """Interactively merge references"""
    try:
        virtool_cli.merge.run(
            source_src_path,
            target_src_path,
            resume,
            in_place,
            output,
            output_unknown_isolates,
            unknown_output,
        )
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


if __name__ == "__main__":
    ref()
