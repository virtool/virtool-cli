from pathlib import Path
import click

import virtool_cli.build
import virtool_cli.divide
import virtool_cli.init
import virtool_cli.isolate
import virtool_cli.repair
import virtool_cli.taxid


@click.group("ref")
def ref():
    """Commands related to reference databases"""
    pass


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a database reference directory",
)
@click.option(
    "-o",
    "--output",
    default="reference.json",
    help="the output path for the reference.json file, defaults to present working directory",
)
@click.option(
    "-V",
    "--version",
    default=None,
    type=str,
    help="the version string to include in the reference.json file",
)
@click.option(
    "-i", "--indent", 
    is_flag=True,
    help="prints the reference.json file with whitespace, for readability")
@click.option(
    "-v", "--verbose", 
    is_flag=True,
    help="prints additional debug information to the console")
def build(src_path, output, indent, version, verbose):
    """Build a Virtool reference JSON file from a reference directory."""
    try:
        virtool_cli.build.run(Path(src_path), Path(output), indent, version)
    except FileNotFoundError as e:
        message = f"Could not find reference directory at '{src_path}'"
        if verbose:
            message += f"\n   {e}"
        click.echo(
            click.style("ERROR: ", fg='red') + message
        )
    except NotADirectoryError as e:
        message = f"Reference directory at '{src_path}' is not formatted correctly \n"
        message += f"  TIP: Try using the '--verbose' flag to identify directory content problems \n"
        if verbose:
            message += f"   {e}"
        click.echo(
            click.style("ERROR: ", fg='red') + message
        )


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to an input 'reference.json' file",
)
@click.option(
    "-o",
    "--output",
    default="src",
    type=str,
    help="the output path for a divided reference directory tree",
)
@click.option(
    "-v", "--verbose", 
    is_flag=True,
    help="prints additional debug information to the console")
def divide(src_path, output, verbose):
    """Divide a reference.json file from Virtool into a reference directory tree."""
    try:
        if not src_path.endswith(".json"):
            raise TypeError
        virtool_cli.divide.run(Path(src_path), Path(output))
    except TypeError as e:
        message = f"'{src_path}' is not a valid JSON file."
        click.echo(
            click.style("  ERROR: ", fg='red') + message)
    except FileNotFoundError as e:
        message = f"'{src_path}' is not a valid file path."
        if verbose:
            message += f"\n  {e}"
        click.echo(
            click.style("  ERROR: ", fg='red') + message)

@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to the input reference.json file",
)
@click.option(
    "-f", "--force_update", 
    is_flag=True,
    help="retrieves taxid from bank regardless of extant information")
@click.option(
    "-v", "--verbose", 
    is_flag=True,
    help="prints additional debug information to the console")
def taxid(src_path, force_update, verbose):
    """Fetch and write taxid for all OTU in given reference directory."""
    try:
        virtool_cli.taxid.run(Path(src_path), force_update)
    except Exception as e:
        message = "Reference directory is invalid"
        if verbose:
            message += f":\n   {e}"
        click.echo(
            click.style("ERROR", fg='red') + message
        )


@ref.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
@click.option(
    "-v", "--verbose", 
    is_flag=True,
    help="prints additional debug information to the console")
def isolate(src_path, verbose):
    """Fetch new isolates for all OTU in a given reference directory."""
    try:
        path = Path(src_path)
        if not path.exists():
            raise FileNotFoundError
        if not path.is_dir():
            raise NotADirectoryError
        virtool_cli.isolate.run(path, verbose=verbose)
    except FileNotFoundError as e:
        message = f"'{src_path}' is not a valid file path."
        if verbose:
            message += f":\n   {e}"
        click.echo(
            click.style("ERROR: ", fg='red') + message
        )
    except NotADirectoryError as e:
        message = f"'{src_path}' is not a valid reference directory \n"
        if verbose:
            message += f"\n   {e}"
        click.echo(
            click.style("ERROR: ", fg='red') + message
        )


@ref.command()
@click.option(
    "-src", "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
def repair(src_path):
    """Fix every OTU in a given reference directory."""
    try:
        virtool_cli.repair.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo(f"{src_path} is not a valid reference directory")

@ref.command()
@click.option(
    "-repo", "--repo_path",
    required=True,
    type=str,
    help="path to a reference repository directory, one level up from source",
)
@click.option(
    "-url", "--api_url",
    type=str,
    required=True,
    envvar='CALLER_API_URL',
    help="URL leading to a tarfile containing caller workflows",
)
@click.option(
    "-cache", "--cache_path",
    envvar='CACHE_PATH',
    required=False,
    type=str,
    help="path to a cache directory, allows for runtime customization",
)
@click.option(
    "-gh", "--github_token",
    type=str,
    default='',
    required=False,
    envvar='GH_PAT',
    help="GitHub Personal Access Token, for more API access permissions as needed",
)
def init(repo_path: str, api_url: str, cache_path:str='', github_token=None):
    """
    Initialize an empty Virtool reference repo for uploading to Github, 
    including workflow callers
    """
    if not cache_path:
        cache_path = Path.home() / '.virtool'
    else:
        cache_path = Path(cache_path)
    if not cache_path.exists():
        cache_path.mkdir(exist_ok=True)
    
    virtool_cli.init.run(
        repo_path=Path(repo_path), 
        api_url=api_url,
        cache_path=cache_path, 
        gh_token=github_token)


if __name__ == "__main__":
    ref()
