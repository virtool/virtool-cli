import click
from virtool_cli.build import build as run_build
from virtool_cli.divide import divide as run_divide
from virtool_cli.taxid import taxid as run_taxid


@click.group()
def cli():
    pass


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to the database src directory")
@click.option("-o", "--output", default="reference.json", help="the output path for the reference.json file")
@click.option("-i", "--indent", is_flag=True)
@click.option("-V", "--version", default=None, type=str,
              help="the version string to include in the reference.json file")
def build(src_path, output, indent, version):
    """
    Build a Virtool reference JSON file from a data directory.

    """
    try:
        run_build(src_path, output, indent, version)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid src directory")


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to the input reference.json file")
@click.option("-o", "--output", default="src", type=str, help="the output path for divided source directory tree")
def divide(src_path, output):
    """
    Divide a reference.json file from Virtool into a src tree.

    """
    try:
        if not src_path.endswith(".json"):
            raise TypeError
        run_divide(src_path, output)
    except (TypeError, FileNotFoundError):
        click.echo("Specified reference file either does not exist or is not a proper JSON file")


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to the input reference.json file")
@click.option("-f", "--force_update", is_flag=True)
def taxid(src_path, force_update):
    """
    Fetch taxid for all OTU in given src directory

    """
    try:
        run_taxid(src_path, force_update)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid src directory")


if __name__ == "__main__":
    cli()
