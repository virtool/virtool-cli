import click

from virtool_cli.ref import ref
from virtool_cli.hmm import hmm
from virtool_cli.catalog.cli import catalog


@click.group()
def entry():
    """
    ex. virtool ref --help

    note: virtool hmm is has additional bioconda dependencies,
    see documentation for installation instructions
    """
    pass


entry.add_command(ref)
entry.add_command(hmm)
entry.add_command(catalog)

entry()
