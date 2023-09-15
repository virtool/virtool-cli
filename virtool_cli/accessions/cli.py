from pathlib import Path
import click

from virtool_cli.accessions.checkup import run as run_checkup
from virtool_cli.accessions.initialize import run as run_init
from virtool_cli.accessions.update import run as run_update
from virtool_cli.accessions.repair import run as run_repair
from virtool_cli.accessions.exclude import run as run_exclude

ERROR_MESSAGE = click.style("ERROR: ", fg='red')


@click.group("acc")
def acc():
    """
    Commands related to accession catalogues.
    """
    pass

@acc.command()
@click.option(
    "-src",
    "--src_path",
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
@click.option('--debug/--no-debug', default=False)
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
            src=Path(src_path),
            catalog=Path(catalog_path),
            debugging=debug
        )
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + f"{src_path} is not a valid reference directory")


@acc.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=str,
    default='.cache/catalog',
    help="the path to a catalog directory",
)
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option('--debug/--no-debug', default=False)
def update(src_path, catalog_path, debug):
    """Generate a catalog of all included accessions in a src directory"""    
    src_dir = Path(src_path) 
    catalog_dir = Path(catalog_path)

    if not src_dir.exists():
        click.echo(ERROR_MESSAGE + f"Source directory does not exist at {src_path}")
        return
    
    if not catalog_dir.exists():
        click.echo(ERROR_MESSAGE + f"Catalog directory does not exist at {catalog_path}")
        return

    try:
        run_update(
            src=src_dir,
            catalog=catalog_dir,
            debugging=debug
        )
    except (FileNotFoundError, NotADirectoryError):
        click.echo(
            ERROR_MESSAGE + \
            f"{src_path} is not a valid reference directory" + \
            f"OR {catalog_path} is not a valid catalog directory")

@acc.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=str,
    help="the path to a catalog directory",
)
@click.option('--debug/--no-debug', default=False)
def checkup(catalog_path, debug):
    """Run a check on the accession catalogue for outstanding issues."""

    if not Path(catalog_path).exists():
        click.echo(ERROR_MESSAGE + f"Catalog directory does not exist at {catalog_path}")

    try:
        run_checkup(Path(catalog_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(ERROR_MESSAGE + f"{catalog_path} is not a valid catalog directory")

@acc.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=str,
    help="the path to a catalog directory",
)
@click.option('--debug/--no-debug', default=False)
def repair(catalog_path, debug):
    """Run a check on the accession catalogue for outstanding issues."""

    if not Path(catalog_path).exists():
        click.echo(ERROR_MESSAGE + f"Catalog directory does not exist at {catalog_path}")

    try:
        run_repair(Path(catalog_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(ERROR_MESSAGE + f"{catalog_path} is not a valid catalog directory")

@acc.command()
@click.option(
    "-cat",
    "--catalog_path",
    required=True,
    type=str,
    help="the path to a catalog directory",
)
@click.option('--debug/--no-debug', default=False)
def exclude(catalog_path, debug):
    """Run a check on the accession catalogue for outstanding issues."""

    if not Path(catalog_path).exists():
        click.echo(ERROR_MESSAGE + f"Catalog directory does not exist at {catalog_path}")

    try:
        run_exclude(Path(catalog_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError):
        click.echo(ERROR_MESSAGE + f"{catalog_path} is not a valid catalog directory")