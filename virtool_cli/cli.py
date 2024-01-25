import click

from virtool_cli.hmm import hmm
from virtool_cli.ref.cli import ref


@click.group()
def entry():
    """ex. virtool ref --help

    note: virtool hmm is has additional bioconda dependencies,
    see documentation for installation instructions
    """


entry.add_command(ref)
entry.add_command(hmm)

entry()
