import click

import virtool_cli.build


@click.command()
@click.option('--src-path', default=1, help='Number of greetings)
def build(src_path):
    """
    Build a Virtool reference JSON file from a data directory.

    """
    virtool_cli.build.build(src_path)
