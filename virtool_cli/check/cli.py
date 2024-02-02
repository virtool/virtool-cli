from pathlib import Path

import click

from virtool_cli.check.otu import check_otu
from virtool_cli.check.reference import check_reference
from virtool_cli.options import debug_option, path_option
from virtool_cli.utils.logging import configure_logger, error_message_style


@click.group("check")
def check():
    """Check the validity of reference data"""


@check.command()
@debug_option
@path_option
def reference(path: Path, debug: bool):
    """Checks the validity of a reference source"""
    configure_logger(debug)

    try:
        check_reference(path)
    except FileNotFoundError:
        click.echo(f"{error_message_style}{path!s} could not be found", err=True)
    except NotADirectoryError:
        click.echo(f"{error_message_style}{path!s} is not a directory", err=True)


@check.command()
@click.argument(
    "otu_path",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@debug_option
def otu(otu_path: Path, debug: bool):
    """Checks the validity of a reference source"""
    configure_logger(debug)

    try:
        check_otu(otu_path)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            f"{error_message_style}{otu_path!s} is not a valid OTU directory",
            err=True,
        )
