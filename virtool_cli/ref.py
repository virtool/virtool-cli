import sys
from pathlib import Path
import click
import structlog

import virtool_cli.build
import virtool_cli.divide
import virtool_cli.isolate
import virtool_cli.repair
import virtool_cli.taxid
from virtool_cli.migrate import run as run_migrate
from virtool_cli.update import run as run_update
from virtool_cli.accessions.catalog import run as run_catalog

ERROR_MESSAGE = click.style("ERROR: ", fg='red')

logger = structlog.get_logger()
    
@click.group("ref")
def ref():
    """
    Commands related to reference files.
    """
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
    type=str,
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
    """Build a Virtool reference JSON file from a source directory."""
    build_logger = logger.bind(command='build', src_path=src_path, out_path=output)

    if not Path(src_path).exists():
        build_logger.error('Directory not found at src_path')
        return
    
    try:
        virtool_cli.build.run(Path(src_path), Path(output), indent, version)
    except (FileNotFoundError, NotADirectoryError):
        build_logger.exception('Source directory has critical errors')
    except:
        build_logger.exception('Unexpected exception')


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
    except (TypeError, FileNotFoundError) as e:
        click.echo(
            ERROR_MESSAGE + "Specified reference file either does not exist or is not a proper JSON file"
        )
        click.echo(e)


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
        click.echo("Scanning OTUs for taxon IDs...")
        virtool_cli.taxid.run(Path(src_path), force_update)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + "Not a valid reference directory")
        click.echo(e)


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
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


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
    if not Path(src_path).exists():
        logger.critical('Source directory does not exist')

    try:
        virtool_cli.repair.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + "Not a valid reference directory")
        click.echo(e)


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option('--debug/--no-debug', default=False)
def migrate(src_path, debug):
    """Convert a reference directory from v1.x to v2.x"""
    if not Path(src_path).exists():
        logger.critical('Source directory does not exist')

    try:
        run_migrate(Path(src_path))
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + "Not a valid reference directory")
        click.echo(e)


@ref.command()
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
def catalog(src_path, catalog_path, debug):
    """Update or generate a catalog of all included accessions in a src directory"""
    if not Path(src_path).exists():
        click.echo(ERROR_MESSAGE + 'Source directory does not exist')
        return
    
    catalog_dir = Path(catalog_path)
    if not catalog_dir.exists():
        catalog_dir.mkdir(parents=True)

    try:
        run_catalog(
            src_path=Path(src_path),
            catalog_path=Path(catalog_path),
            debugging=debug
        )
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + "Not a valid reference directory")
        logger.exception(e)


@ref.command()
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
@click.option('--debug/--no-debug', default=False)
def update(src_path, catalog_path, debug):
    """Fetch new isolates for all OTU in a given reference directory."""
    if not Path(catalog_path).exists():
        click.echo("Not a valid catalog directory")
        return

    try:
        run_update(Path(src_path), Path(catalog_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)

if __name__ == "__main__":
    ref()
