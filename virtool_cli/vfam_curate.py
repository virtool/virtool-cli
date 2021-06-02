import os
import sys


from Bio import SeqIO
from pathlib import Path


SEQUENCE_MIN_LENGTH = 1


def get_input_paths(src_path: Path) -> list:
    """
    Takes in a path to directory containing input fasta files, returns a list of paths to the fasta files

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

    :return: no_phages, list of filtered records
    """
    no_phages = list()

    for input_path in input_paths:
        with open(input_path) as handle:

            for record in SeqIO.parse(handle, "fasta"):
                if "phage" not in record.description:
                    no_phages.append(record)

    return no_phages


def remove_dupes(no_phages: list, output: Path) -> Path:
    """
    Removes duplicates in no_phages list, writes all records in list to output

    verifies length of each sequence is longer than SEQUENCE_MIN_LENGTH

    :param no_phages: list of records from all protein files without keyword "phage"
    :param output: Path to output directory for profile HMMs and intermediate files
    :return: Path to filtered file
    """
    record_ids = []
    records_to_output = []

    for record in no_phages:
        if record.id not in record_ids and len(record.seq) > SEQUENCE_MIN_LENGTH:
            records_to_output.append(record)
            record_ids.append(record.id)

    SeqIO.write(records_to_output, output, "fasta")

    return output

