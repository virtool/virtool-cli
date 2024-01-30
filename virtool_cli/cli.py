import click

from virtool_cli.ref.cli import ref


@click.group()
def entry():
    """Manage Virtool datasets"""


entry.add_command(ref)

entry()
