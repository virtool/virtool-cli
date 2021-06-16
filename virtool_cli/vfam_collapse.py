import subprocess

from pathlib import Path
from Bio import SeqIO


def generate_clusters(curated_fasta: Path, prefix, fraction_cov, fraction_id: float) -> Path:
    """
    Takes in a fasta file, minimum fraction coverage, minimum fraction identity, calls cd-hit to cluster data

    cd-hit collapses the input sequences into non-redundant representatives at the specified levels

    :param curated_fasta: fasta file from vfam_curate step, to have cd-hit performed on it
    :param prefix: Prefix for intermediate and result files
    :param fraction_cov: Fraction coverage for cd-hit step
    :param fraction_id: Fraction ID for cd-hit step
    :return: output_path, path to a file containing cluster information created at cd-hit step
    """
    output_name = "clustered_records.faa"
    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = curated_fasta.parent / output_name

    cd_hit_cmd = ["cd-hit", "-i", str(curated_fasta), "-o", str(output_path), "-s", str(fraction_id)]

    if fraction_cov:
        cd_hit_cmd += ["-c", str(fraction_cov)]

    subprocess.run(cd_hit_cmd, check=True, shell=False)

    return output_path


def polyprotein_name_check(clustered_records: Path, prefix) -> Path:
    """
    Removes any record with "polyprotein" found in its description

    :param clustered_records: file containing clustered fasta information from generate clusters step
    :param prefix: Prefix for intermediate and result files
    :return: output_path, a file containing all clustered fasta records with polyproteins removed
    """
    output_name = "no_named_polyproteins.faa"
    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = clustered_records.parent / Path(output_name)

    filtered_records = []

    for record in SeqIO.parse(clustered_records, "fasta"):
        if "polyprotein" not in record.description:
            filtered_records.append(record)

    SeqIO.write(filtered_records, output_path, "fasta")

    return output_path


def all_by_all_blast(clustered_fasta: Path, prefix, num_cores: int) -> Path:
    """
    Takes a clustered fasta file as input, and formats file to be a BLAST protein database

    runs BLAST on the file to itself as the BLAST-formatted database

    :param clustered_fasta: fasta file from previous cd-hit step
    :param prefix: Prefix for intermediate and result files
    :param num_cores: number of cores to use at all_by_all blast step
    :return: tab-delimited formatted BLAST results file
    """
    format_db_cmd = ["makeblastdb", "-in", clustered_fasta, "-dbtype", "prot"]
    subprocess.run(format_db_cmd, check=True, shell=False)

    blast_name = "blast.br"
    if prefix:
        blast_name = f"{prefix}_{blast_name}"

    blast_results = clustered_fasta.parent / Path(blast_name)

    blast_cmd = ["blastp", "-query", clustered_fasta, "-out", blast_results, "-db", clustered_fasta,
                 "-outfmt", "6", "-num_threads", str(num_cores)]

    subprocess.run(blast_cmd, check=True, shell=False)

    return blast_results
