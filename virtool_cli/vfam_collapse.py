import subprocess
from pathlib import Path


NUM_CORES = 8


def generate_clusters(curated_fasta: Path, output: Path, fraction_cov: float, fraction_id: float) -> Path:
    """
    Takes in a fasta file, minimum fraction coverage, minimum fraction identity, calls cd-hit to cluster data

    cd-hit collapses the input sequences into non-redundant representatives at the specified levels

    :param curated_fasta: fasta file from vfam_curate step, to have cd-hit performed on it
    :param output: Path to output directory
    :param fraction_cov: Fraction coverage for cd-hit step
    :param fraction_id: Fraction ID for cd-hit step
    :return: cluster_file, a file containing cluster information created at cd-hit step
    """
    output_path = output / Path("clustered_records")

    if not fraction_cov:
        coverage_args = ""
    else:
        coverage_args = "-c %s" % fraction_cov

    cd_hit_cmd = "cd-hit -i %s -o %s %s -s %s" % (
        curated_fasta, output_path, coverage_args, fraction_id)
    subprocess.run(cd_hit_cmd.split())

    cluster_file = str(output_path) + ".clstr"

    return Path(cluster_file)


def all_by_all_blast(clustered_fasta_file) -> Path:
    """
    Takes a clustered fasta file as input, and formats file to be a BLAST protein database

    runs BLAST on the file to itself as the BLAST-formatted database

    :param clustered_fasta_file:
    :return: tab-delimited formatted BLAST results file
    """

    format_database_cmd = "makeblastdb -in %s -dbtype prot" % clustered_fasta_file
    subprocess.run(format_database_cmd.split())

    blast_results_file = Path(clustered_fasta_file).parent / "all_by_all_blast_results.br"
    all_by_all_blast_cmd = "blastp -query %s -out %s -db %s -outfmt 6 -num_threads %s" % (
        clustered_fasta_file, blast_results_file, clustered_fasta_file, NUM_CORES)
    subprocess.run(all_by_all_blast_cmd.split())

    return blast_results_file