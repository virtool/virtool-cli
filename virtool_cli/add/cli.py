from pathlib import Path
import click
import logging

from virtool_cli.add.accession import run as run_add_accession

# from virtool_cli.add.otu import run as run_add_otu


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
    help="An accession from the NCBI Taxonomy database",
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


#
#
# @add.command()
# @click.option(
#     "-otu",
#     "--otu_path",
#     required=True,
#     type=click.Path(exists=True, file_okay=False, path_type=Path),
#     help="the path to a single OTU directory",
# )
# @click.option(
#     "-cat",
#     "--catalog_path",
#     required=True,
#     type=click.Path(exists=True, file_okay=False, path_type=Path),
#     default=".cache/catalog",
#     help="the path to a catalog directory",
# )
# @click.option("--evaluate/--no-evaluate", default=False, help="Enable auto-filtering")
# @click.option("--debug/--no-debug", default=False, help="Enable debugging logs")
# def otu(otu_path, catalog_path, evaluate, debug):
#     """Fetch new sequences and isolates for a given OTU directory."""
#
#     try:
#         run_update_single(
#             otu_path, catalog_path, auto_evaluate=evaluate, debugging=debug
#         )
#
#     except (FileNotFoundError, NotADirectoryError) as e:
#         click.echo("Not a valid reference directory")
#         click.echo(e)


if __name__ == "__main__":
    update()
