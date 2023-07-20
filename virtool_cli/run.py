import click

from virtool_cli.ref import ref
from virtool_cli.hmm import hmm


class OrderCommands(click.Group):
  """
  Modifies Group to list commands in the order of definition, 
  instead of alphanumeric order
  """
  def list_commands(self, ctx: click.Context) -> list[str]:
    return list(self.commands)

@click.group(cls=OrderCommands)
def cli():
    """
    ex. virtool ref --help

    note: virtool hmm is has additional bioconda dependencies,
    see documentation for installation instructions 
    """
    pass


cli.add_command(ref)
cli.add_command(hmm)
cli()