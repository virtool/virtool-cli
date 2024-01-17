import subprocess
import sys
import click
from pathlib import Path

from virtool_cli.vfam.cmd import run

ERROR_MESSAGE = click.style("ERROR: ", fg="red")


@click.group("hmm")
def hmm():
    """
    Commands related to Hidden Markov Models.

    Requires bioconda packages: cd-hit, hmmer, blast, mcl, muscle
    """
    pass


@hmm.command()
@click.option(
    "--src-path",
    default=None,
    help="Path to input reference directory if not gathering from NCBI",
    type=click.Path(exists=True)
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
    src_path: Path,
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
    """
    Build profile HMMs from FASTAs
    """
    try:
        check_vfam_dependencies()
    except (FileNotFoundError, PermissionError) as e:
        click.echo(ERROR_MESSAGE + "Missing external program dependency")
        click.echo(e)
        sys.exit(1)

    try:
        run(
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
        click.echo(ERROR_MESSAGE + "Not a valid reference directory.")


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


if __name__ == "__main__":
    hmm()
