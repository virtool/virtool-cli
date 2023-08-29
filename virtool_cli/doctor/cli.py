from pathlib import Path
import click
# import structlog

from virtool_cli.doctor.fix_reference import run as run_fix_all
# from virtool_cli.doctor.update import run as run_update_single

@click.group("doc")
def doc():
    """
    Commands related to reference repairs
    """
    pass

@doc.command()
@click.option(
    "-src", "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
# @click.option(
#     "-cat",
#     "--catalog_path",
#     required=True,
#     type=str,
#     default='.cache/catalog',
#     help="the path to a catalog directory",
# )
@click.option('--debug/--no-debug', default=False)
def reference(src_path, debug):
    """Fetch new sequences and isolates for all OTU in a given reference directory."""
    # if not Path(catalog_path).exists():
    #     click.echo("Not a valid catalog directory")
    #     return

    try:
        run_fix_all(Path(src_path), debugging=debug)
    except (FileNotFoundError, NotADirectoryError) as e:
        click.echo("Not a valid reference directory")
        click.echo(e)

# @doc.command()
# @click.option(
#     "-otu", "--otu_path",
#     required=True,
#     type=str,
#     help="the path to a single OTU directory",
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
# def otu(otu_path, catalog_path, debug):
#     """Fetch new sequences and isolates for a given OTU directory."""
#     if not Path(catalog_path).exists():
#         click.echo("Not a valid catalog directory")
#         return

#     try:
#         run_update_single(Path(otu_path), Path(catalog_path), debugging=debug)
#     except (FileNotFoundError, NotADirectoryError) as e:
#         click.echo("Not a valid reference directory")
#         click.echo(e)


if __name__ == "__main__":
    doc()