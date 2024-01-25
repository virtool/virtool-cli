from pathlib import Path

import click

debug_option = click.option("--debug", is_flag=True, help="Show debug logs")
"""A click option for enabling debug logs."""

path_option = click.option(
    "--path",
    default=".",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to the reference repository",
)
"""A click option for the path to the reference repository."""
