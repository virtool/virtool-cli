import os
import sys


from Bio import SeqIO
from pathlib import Path


def get_input_paths(src_path: Path) -> list:
    """
    Takes in a path to directory containing input fasta files, returns a list of paths to the fasta files

    :param src_path: Path to input source directory containing unfiltered fasta files
    :return: input_paths, list of paths to input protein files if any files are found
    """
    input_files = os.listdir(src_path)
    input_paths = [src_path / file for file in input_files]

    if input_paths:
        return input_paths
    else:
        print("Error: no files found in input directory")
        sys.exit(1)


def remove_phages(input_paths: list) -> list:
    """
    Parses through records in input files, appends records to filtered_sequences if keyword "phage" not found in record

    :param input_paths: list of paths to input fasta files
    :return: no_phages, list of curated records with all phage records removed
    """
    no_phages = list()

    for input_path in input_paths:
        with open(input_path) as handle:

            for record in SeqIO.parse(handle, "fasta"):
                if "phage" not in record.description:
                    no_phages.append(record)

    return no_phages


def group_input_paths(input_paths):
    """
    Takes input paths as input and created list of all records found in each path

    :param input_paths: list of paths to input files
    :return: list of all records found in input paths
    """
    records = []

    for input_path in input_paths:
        with open(input_path) as handle:

            for record in SeqIO.parse(handle, "fasta"):
                records.append(record)

    return records


def remove_dupes(no_phages: list, output: Path, sequence_min_length: int) -> Path:
    """
    Removes duplicates in no_phages list, writes all records in list to output

    verifies length of each sequence is longer than SEQUENCE_MIN_LENGTH

    :param no_phages: list of records from all protein files without keyword "phage"
    :param output: Path to output directory for profile HMMs and intermediate files
    :param sequence_min_length: Minimum sequence length for a record to be included in the input
    :return: Path to curated fasta file without repeats or phages
    """
    record_seqs = []
    record_ids = []

    records_to_output = []

    for record in no_phages:
        if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
            if record.description:
                records_to_output.append(record)
                record_seqs.append(record.seq)
                record_ids.append(record.id)

    output_dir = output / Path("cluster_files")

    if not output_dir.exists():
        output_dir.mkdir()

    output_path = output_dir / Path("curated_records.faa")

    SeqIO.write(records_to_output, Path(output_path), "fasta")

    return output_path



