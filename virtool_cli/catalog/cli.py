from pathlib import Path
import click

from virtool_cli.catalog.checkup import run as run_checkup
from virtool_cli.catalog.initialize import run as run_init
from virtool_cli.catalog.update import run as run_update
from virtool_cli.catalog.repair import run as run_repair
from virtool_cli.catalog.exclude import run as run_exclude

ERROR_MESSAGE = click.style("ERROR: ", fg="red")


@click.group("catalog")
def catalog():
    """
    Commands related to accession catalogues.
    """
    pass


@catalog.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    default=".cache/catalog",
    help="the path to a catalog directory",
)
@click.option("--debug/--no-debug", default=False)
def init(src_path, catalog_path, debug):
    """Generate a catalog of all included accessions in a src directory"""

    if not Path(src_path).exists():
        click.echo(ERROR_MESSAGE + f"Source directory does not exist at {src_path}")
        return

    catalog_dir = Path(catalog_path)
    if not catalog_dir.exists():
        catalog_dir.mkdir(parents=True)

    try:
        run_init(
            src_path=Path(src_path), catalog_path=Path(catalog_path), debugging=debug
        )
    except (FileNotFoundError, NotADirectoryError):
        click.echo(ERROR_MESSAGE + f"{src_path} is not a valid reference directory")


@catalog.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=".cache/catalog",
    help="the path to a catalog directory",
)
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option("--debug/--no-debug", default=False)
def update(src_path, catalog_path, debug):
    """Generate a catalog of all included accessions in a src directory"""

    try:
        run_update(src_path=src_path, catalog_path=catalog_path, debugging=debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            ERROR_MESSAGE
            + f"{str(src_path)} is not a valid reference directory"
            + f"OR {str(catalog_path)} is not a valid catalog directory"
        )


@catalog.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a catalog directory",
)
@click.option("--debug/--no-debug", default=False)
def checkup(catalog_path, debug):
    """Run a check on the accession catalogue for outstanding issues."""

    try:
        run_checkup(Path(catalog_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(ERROR_MESSAGE + f"{catalog_path} is not a valid catalog directory")


@catalog.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a catalog directory",
)
@click.option("--debug/--no-debug", default=False)
def repair(catalog_path, debug):
    """Run a check on the accession catalogue for outstanding issues."""

    try:
        run_repair(Path(catalog_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(ERROR_MESSAGE + f"{catalog_path} is not a valid catalog directory")


@catalog.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a catalog directory",
)
@click.option("--debug/--no-debug", default=False)
def exclude(catalog_path, debug):
    """Run a check on the accession catalogue for outstanding issues."""

    try:
        run_exclude(catalog_path, debugging=debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            ERROR_MESSAGE + f"{str(catalog_path)} is not a valid catalog directory"
        )
