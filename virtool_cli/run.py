import subprocess
import sys
import click
import virtool_cli.vfam
import virtool_cli.build
import virtool_cli.divide
import virtool_cli.isolate
import virtool_cli.repair
import virtool_cli.taxid
import virtool_cli.merge_refs

from pathlib import Path
from virtool_cli.vfam_console import console


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "-src",
    "--src-path",
    default=None,
    help="Path to input reference directory if not gathering from NCBI",
)
@click.option(
    "-o",
    "--output",
    required=True,
    help="Path to output directory for profile HMMs and " "intermediate files",
)
@click.option(
    "-p", "--prefix", default=None, help="Prefix for intermediate and result files"
)
@click.option(
    "--sequence-min-length",
    default=1,
    help="Minimum sequence length to be included in input",
)
@click.option(
    "--no-named-phages",
    default=False,
    help="Filter out phage sequences based on record description",
)
@click.option(
    "--fraction-coverage", default=None, help="Fraction coverage for cd-hit step"
)
@click.option("--fraction-id", default=1.0, help="Fraction ID for cd-hit step")
@click.option(
    "-n",
    "--num-cores",
    default=8,
    help="Number of cores to be used in all by all blast step",
)
@click.option(
    "--no-named-polyproteins",
    default=False,
    help="Filter out polyprotein sequences based on " "record description",
)
@click.option(
    "--inflation-num", default=None, help="Inflation number to be used in mcl step"
)
@click.option(
    "--filter-clusters", default=False, help="Filter clustered fasta files on coverage"
)
@click.option(
    "-min-seqs",
    "--min-sequences",
    default=2,
    help="Filter out clusters with fewer records than " "min-sequences",
)
def vfam(
    src_path: str,
    output: str,
    prefix: str,
    sequence_min_length: int,
    no_named_phages: bool,
    fraction_coverage: float,
    fraction_id: float,
    num_cores: int,
    no_named_polyproteins: bool,
    inflation_num: int,
    filter_clusters: bool,
    min_sequences: int,
):
    """Build profile HMMS from fasta."""
    try:
        check_vfam_dependencies()
    except (FileNotFoundError, PermissionError):
        console.print("Missing external program dependency.", style="red")
        sys.exit(1)

    try:
        virtool_cli.vfam.run(
            src_path,
            Path(output),
            prefix,
            sequence_min_length,
            no_named_phages,
            fraction_coverage,
            fraction_id,
            num_cores,
            no_named_polyproteins,
            inflation_num,
            filter_clusters,
            min_sequences,
        )
    except (FileNotFoundError, NotADirectoryError):
        console.print("Not a valid reference directory.", style="red")


def check_vfam_dependencies():
    """Check external dependencies for VFam pipeline"""
    subprocess.run(["hmmstat", "-h"])
    subprocess.run(["cd-hit", "-h"])
    subprocess.run(["makeblastdb", "-h"])
    subprocess.run(["blastp", "-h"])
    subprocess.run(["mcxload", "-h"])
    subprocess.run(["mcl", "-h"])
    subprocess.run(["muscle", "-h"])
    subprocess.run(["hmmbuild", "-h"])


@cli.command()
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
    help="the output path for a reference.json file",
)
@click.option("-i", "--indent", is_flag=True)
@click.option(
    "-V",
    "--version",
    default=None,
    type=str,
    help="the version string to include in the reference.json file",
)
def build(src_path, output, indent, version):
    """Build a Virtool reference JSON file from a reference directory."""
    try:
        virtool_cli.build.run(Path(src_path), Path(output), indent, version)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a input reference.json file",
)
@click.option(
    "-o",
    "--output",
    default="src",
    type=str,
    help="the output path for a divided reference directory tree",
)
def divide(src_path, output):
    """Divide a reference.json file from Virtool into a reference directory tree."""
    try:
        if not src_path.endswith(".json"):
            raise TypeError
        virtool_cli.divide.run(Path(src_path), Path(output))
    except (TypeError, FileNotFoundError):
        click.echo(
            "Specified reference file either does not exist or is not a proper JSON file"
        )


@cli.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to the input reference.json file",
)
@click.option("-f", "--force_update", is_flag=True)
def taxid(src_path, force_update):
    """Fetch taxid for all OTU in given reference directory."""
    try:
        virtool_cli.taxid.run(Path(src_path), force_update)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
def isolate(src_path):
    """Fetch new isolates for all OTU in a given reference directory."""
    try:
        virtool_cli.isolate.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option(
    "-src",
    "--src_path",
    required=True,
    type=str,
    help="the path to a reference directory",
)
def repair(src_path):
    """Fix every OTU in a given reference directory."""
    try:
        virtool_cli.repair.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option(
    "-s",
    "--source-src-path",
    required=True,
    type=str,
    help="Path to src for reference to copy isolates from",
)
@click.option(
    "-t",
    "--target-src-path",
    required=True,
    type=str,
    help="Path to src for reference to write to",
)
@click.option(
    "-r", "--resume", is_flag=True, help="Resume in-progress merge using cache"
)
@click.option(
    "-i",
    "--in-place",
    is_flag=True,
    help="Place new isolates directly into target reference",
)
@click.option(
    "-o", "--output", default="merged_ref", help="Name for merged reference directory"
)
def merge_refs(source_src_path, target_src_path, resume, in_place, output):
    try:
        virtool_cli.merge_refs.run(
            source_src_path, target_src_path, resume, in_place, output
        )
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


if __name__ == "__main__":
    cli()
