import click
from virtool_cli.ref import ref
from virtool_cli.hmm import hmm


@click.group()
def cli():
    pass


cli.add_command(ref)
cli.add_command(hmm)
cli()
