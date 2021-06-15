from pathlib import Path

import click

import virtool_cli.vfam
import virtool_cli.build
import virtool_cli.divide
import virtool_cli.isolate
import virtool_cli.repair
import virtool_cli.taxid


@click.group()
def cli():
    pass


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="Path to input reference directory")
@click.option("-o", "--output", required=True, help="Path to output directory for profile HMMs and intermediate files")
@click.option("-p", "--prefix", default=None, help="Prefix for intermediate and result files")
@click.option("-seq_min", "--sequence_min_length", default=1, help="Minimum sequence length to be included in "
                                                                   "input")
@click.option("-phage", "phage_name_check", default=False,
              help="Filter out phage sequences based on record description")
@click.option("-f_cov", "--fraction_coverage", default=None, help="Fraction coverage for cd-hit step")
@click.option("-f_id", "--fraction_id", default=1.0, help="Fraction ID for cd-hit step")
@click.option("-cores", "--num_cores", default=8, help="Number of cores to be used in all by all blast step")
@click.option("-polyp", "--polyp_name_check", default=False, help="Filter out polyprotein sequences based on "
                                                                 "record description")
@click.option("-i_num", "--inflation_num", default=None, help="Inflation number to be used in mcl step")
@click.option("-cvg_check", "--filter_on_cvg", default=False, help="Filter clustered fasta files on coverage")
@click.option("-min_seqs", "--min_sequences", default=2, help="Filter out clusters with fewer records than "
                                                              "min_sequences")
def vfam(src_path, output, prefix, sequence_min_length, phage_name_check, fraction_coverage, fraction_id, num_cores,
         polyp_name_check, inflation_num, filter_on_cvg, min_sequences):
    """Build profile HMMS from fasta."""
    try:
        virtool_cli.vfam.run(Path(src_path), Path(output), prefix, sequence_min_length, phage_name_check,
                             fraction_coverage, fraction_id, num_cores, polyp_name_check, inflation_num,
                             filter_on_cvg, min_sequences)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to a database reference directory")
@click.option("-o", "--output", default="reference.json", help="the output path for a reference.json file")
@click.option("-i", "--indent", is_flag=True)
@click.option("-V", "--version", default=None, type=str,
              help="the version string to include in the reference.json file")
def build(src_path, output, indent, version):
    """Build a Virtool reference JSON file from a reference directory."""
    try:
        virtool_cli.build.run(Path(src_path), Path(output), indent, version)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to a input reference.json file")
@click.option("-o", "--output", default="src", type=str, help="the output path for a divided reference directory tree")
def divide(src_path, output):
    """Divide a reference.json file from Virtool into a reference directory tree."""
    try:
        if not src_path.endswith(".json"):
            raise TypeError
        virtool_cli.divide.run(Path(src_path), Path(output))
    except (TypeError, FileNotFoundError):
        click.echo("Specified reference file either does not exist or is not a proper JSON file")


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to the input reference.json file")
@click.option("-f", "--force_update", is_flag=True)
def taxid(src_path, force_update):
    """Fetch taxid for all OTU in given reference directory."""
    try:
        virtool_cli.taxid.run(Path(src_path), force_update)
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to a reference directory")
def isolate(src_path):
    """Fetch new isolates for all OTU in a given reference directory."""
    try:
        virtool_cli.isolate.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


@cli.command()
@click.option("-src", "--src_path", required=True, type=str, help="the path to a reference directory")
def repair(src_path):
    """Fix every OTU in a given reference directory."""
    try:
        virtool_cli.repair.run(Path(src_path))
    except (FileNotFoundError, NotADirectoryError):
        click.echo("Not a valid reference directory")


if __name__ == "__main__":
    cli()
