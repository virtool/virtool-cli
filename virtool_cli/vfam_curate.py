import sys

from Bio import SeqIO
from pathlib import Path
from typing import List, Optional


def get_input_paths(src_path: Path) -> List[Path]:
    """
    Takes in a path to directory containing input fasta files, returns a list of paths to the fasta files

    :param src_path: Path to input source directory containing unfiltered fasta files
    :return: input_paths, list of paths to input files if any files are found
    """
    input_paths = list(src_path.iterdir())

    if input_paths:
        return input_paths
    else:
        print("Error: no files found in input directory")
        sys.exit(1)


def group_input_paths(input_paths: list, no_named_phages: bool) -> list:
    """
    Takes input paths as input and yields records

    filters out records with "phage" in their description if id no_named_phages is True

    :param input_paths: list of paths to input files
    :param no_named_phages: bool that dictates whether phage records are filtered out or not
    :return: list of all records found in input paths
    """
    for input_path in input_paths:
        for record in SeqIO.parse(input_path, "fasta"):
            if no_named_phages:
                if "phage" not in record.description:
                    yield record
            else:
                yield record


def remove_dupes(records: iter, sequence_min_length: int):
    """
    Iterates through records and yields record if sequence isn't in record_seqs and is longer than sequence_min_length

    :param records: iterable of records gathered in group_input_paths()
    :param sequence_min_length: minimum length of sequence to be included in output
    """
    record_seqs = []

    for record in records:
        if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
            record_seqs.append(record.seq)
            yield record


def write_curated_recs(records: iter, output: Path, prefix: Optional[str], sequence_min_length: int) -> Path:
    """
    Removes duplicates in no_phages list, writes all records in list to output

    verifies length of each sequence is longer than SEQUENCE_MIN_LENGTH

    :param records: list of records from all protein files without keyword "phage"
    :param output: Path to output directory for profile HMMs and intermediate files
    :param prefix: Prefix for intermediate and result files
    :param sequence_min_length: Minimum sequence length for a record to be included in the input
    :return: Path to curated fasta file without repeats or phages
    """
    record_seqs = []
    records_to_output = []

    for record in records:
        if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
            if record.description:
                records_to_output.append(record)
                record_seqs.append(record.seq)

    output_dir = output / Path("intermediate_files")

    if not output_dir.exists():
        output_dir.mkdir()

    output_name = "curated_records.faa"
    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = output_dir / Path(output_name)
    SeqIO.write(records_to_output, Path(output_path), "fasta")

    return output_path
