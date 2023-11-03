from pathlib import Path
import click

from virtool_cli.add.accession import run_single as run_add_accession
from virtool_cli.add.accession import run_multiple as run_add_accessions
from virtool_cli.add.otu import run as run_add_otu


@click.group("add")
def add():
    """
    Commands related to new data
    """
    pass


@add.command()
@click.option(
    "-acc",
    "--accession",
    required=True,
    type=str,
    help="An accession from the NCBI Nucleotide database",
)
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
    required=False,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=".cache/catalog",
    help="the path to a catalog directory",
)
@click.option("--debug/--no-debug", default=False, help="Enable debugging logs")
def accession(accession, src_path, catalog_path, debug):
    """Takes an accession and retrieves the data."""
    try:
        run_add_accession(accession, src_path, catalog_path, debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


@add.command()
@click.option(
    "-acc",
    "--accessions",
    required=True,
    type=str,
    help="An accession from the NCBI Nucleotide database",
)
@click.option(
    "-otu",
    "--otu_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a reference directory",
)
@click.option(
    "-cat",
    "--catalog_path",
    required=False,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=".cache/catalog",
    help="the path to a catalog directory",
)
@click.option("--debug/--no-debug", default=False, help="Enable debugging logs")
def accessions(accessions, otu_path, catalog_path, debug):
    """Takes a list of comma-delineated accessions and retrieves the data."""
    try:
        run_add_accessions(accessions, otu_path, catalog_path, debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


@add.command()
@click.option(
    "-taxid",
    "--taxon_id",
    required=True,
    type=int,
    help="A unique identifier from the NCBI Taxonomy database",
)
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
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=".cache/catalog",
    help="the path to a catalog directory",
)
@click.option("--debug/--no-debug", default=False, help="Enable debugging logs")
def otu(taxon_id, src_path, catalog_path, debug):
    """Create a new taxon ID and populate with accessions."""

    try:
        run_add_otu(taxon_id, src_path, catalog_path, debugging=debug)

    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)


if __name__ == "__main__":
    add()
