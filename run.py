import click
import virtool_cli.build


@click.group()
def cli():
    pass


@cli.command()
@click.argument('src_path', nargs=1, required=True)
@click.option("--output", default="reference.json", help="the output path for the reference.json file")
@click.option("-i", "--indent", is_flag=True)
@click.option("-V", "--version", default=None, type=str,
              help="the version string to include in the reference.json file")
def build(src_path, output, indent, version):
    """
    Build a Virtool reference JSON file from a data directory.\n
    SRC_PATH: the path to the database src directory

    """
    virtool_cli.build.build(src_path, output, indent, version)


if __name__ == "__main__":
    cli()
