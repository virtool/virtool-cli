import subprocess
from pathlib import Path
from Bio import SeqIO

def generate_clusters(curated_fasta: Path, output: Path, fraction_cov, fraction_id: float) -> Path:
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

    return output_path


def polyprotein_name_check(clustered_records, output):
    no_polyproteins = []
    polyproteins_removed = output / Path("polyp_names_removed")

    with clustered_records as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if "polyprotein" not in record.description:
                no_polyproteins.append(record)
    with polyproteins_removed as filtered_records:
        SeqIO.write(no_polyproteins, filtered_records, "fasta")

    return polyproteins_removed


def all_by_all_blast(clustered_fasta: Path, output: Path, num_cores) -> Path:
    """
    Takes a clustered fasta file as input, and formats file to be a BLAST protein database

    runs BLAST on the file to itself as the BLAST-formatted database

    :param clustered_fasta: fasta file from previous cd-hit step
    :param output: output directory for project
    :param num_cores: number of cores to use at all_by_all blast step
    :return: tab-delimited formatted BLAST results file
    """
    format_db_cmd = "makeblastdb -in %s -dbtype prot" % clustered_fasta
    subprocess.run(format_db_cmd.split())

    blast_results_file = output / Path("all_by_all.br")

    all_by_all_blast_cmd = "blastp -query %s -out %s -db %s -outfmt 6 -num_threads %s" % (
        clustered_fasta, blast_results_file, clustered_fasta, num_cores)
    subprocess.run(all_by_all_blast_cmd.split())

    return blast_results_file
