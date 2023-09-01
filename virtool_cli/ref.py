import sys
from pathlib import Path
import click
# import structlog

import virtool_cli.build
import virtool_cli.divide
from virtool_cli.doctor.cli import doctor
from virtool_cli.update.cli import update
from virtool_cli.migrate import run as run_migrate
# from virtool_cli.update.update_ref import run as run_update
# from virtool_cli.accessions.catalog import run as run_catalog

ERROR_MESSAGE = click.style("ERROR: ", fg='red')

    
@click.group("ref")
def ref():
    """
    Commands related to reference files.
    """
    pass

ref.add_command(update)
ref.add_command(doctor)

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
@click.option('--debug/--no-debug', default=False)
def build(src_path, output, indent, version, debug):
    """Build a Virtool reference JSON file from a source directory."""
    if not Path(src_path).exists():
        # build_logger.error('Directory not found at src_path')
        click.echo(ERROR_MESSAGE + f'No directory found at {src_path}')
        return
    
    try:
        virtool_cli.build.run(
            Path(src_path), Path(output), indent, version, debug
        )
    except (FileNotFoundError, NotADirectoryError):
        click.echo(ERROR_MESSAGE + f'Source directory has critical errors')
    except:
        click.echo(ERROR_MESSAGE + f'Unexpected exception')


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
@click.option('--debug/--no-debug', default=False)
def divide(src_path, output, debug):
    """Divide a reference.json file from Virtool into a reference directory tree."""
    try:
        if not src_path.endswith(".json"):
            raise TypeError
        virtool_cli.divide.run(
            Path(src_path), Path(output), 
            debug)
    except (TypeError, FileNotFoundError) as e:
        click.echo(
            ERROR_MESSAGE + f"{src_path} either does not exist or is not a proper JSON file"
        )
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
        click.echo(ERROR_MESSAGE + f'No directory found at {src_path}')

    try:
        run_migrate(Path(src_path), debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo(ERROR_MESSAGE + f"{src_path} is not a valid reference directory")
        click.echo(e)


# @ref.command()
# @click.option(
#     "-src",
#     "--src_path",
#     required=True,
#     type=str,
#     help="the path to a reference directory",
# )
# @click.option(
#     "-cat",
#     "--catalog_path",
#     required=True,
#     type=str,
#     default='.cache/catalog',
#     help="the path to a catalog directory",
# )
# @click.option('--debug/--no-debug', default=False)
# def catalog(src_path, catalog_path, debug):
#     """Update or generate a catalog of all included accessions in a src directory"""
#     if not Path(src_path).exists():
#         click.echo(ERROR_MESSAGE + 'Source directory does not exist')
#         return
    
#     catalog_dir = Path(catalog_path)
#     if not catalog_dir.exists():
#         catalog_dir.mkdir(parents=True)

#     try:
#         run_catalog(
#             src_path=Path(src_path),
#             catalog_path=Path(catalog_path),
#             debugging=debug
#         )
#     except (FileNotFoundError, NotADirectoryError) as e:
#         click.echo(ERROR_MESSAGE + "Not a valid reference directory")
#         logger.exception(e)

if __name__ == "__main__":
    ref()
