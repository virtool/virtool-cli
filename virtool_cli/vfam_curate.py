import sys

from Bio import SeqIO
from pathlib import Path
from typing import List
from virtool_cli.vfam_console import console


def get_input_paths(src_path: Path) -> List[Path]:
    """
    Takes in a path to directory containing input FASTA files, returns a list of paths to the FASTA files.

    :param src_path: Path to input source directory containing unfiltered FASTA files
    :return: input_paths, list of paths to input files if any files are found
    """
    input_paths = list(src_path.iterdir())

    if input_paths:
        console.print(f"✔ Retrieved {len(input_paths)} files from input directory.", style="green")
        return input_paths

    console.print("No files found in input directory.", style="red")
    sys.exit(1)


def group_input_paths(input_paths: List[Path], no_named_phages: bool) -> list:
    """
    Takes input_paths from get_input_paths() as input and yields records.

    Filters out records with "phage" in their description if no_named_phages is True.

    :param input_paths: list of paths to input FASTA files
    :param no_named_phages: bool that dictates whether phage records are filtered out by name or not
    :return: list of all records found in input paths
    """
    phage_count = 0
    record_count = 0
    for input_path in input_paths:
        for record in SeqIO.parse(input_path, "fasta"):
            record_count += 1
            if no_named_phages:
                if "phage" in record.description:
                    phage_count += 1
                else:
                    yield record
            else:
                yield record

    console.print(f"✔ Retrieved {record_count} records from {len(input_paths)} input files.", style="green")

    if no_named_phages:
        console.print(f"✔ Filtered out {phage_count} phage records by name.", style="green")


def remove_dupes(records: iter, sequence_min_length: int) -> list:
    """
    Iterates through records and yields record if sequence isn't in record_seqs and is longer than sequence_min_length.

    :param records: iterable of records gathered in group_input_paths()
    :param sequence_min_length: minimum length of sequence to be included in output
    """
    record_seqs = list()
    dupes_count = 0

    for record in records:
        if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
            record_seqs.append(record.seq)
            yield record
        else:
            dupes_count += 1

    console.print(f"✔ Filtered out {dupes_count} duplicate records.", style="green")


def write_curated_recs(records: iter, output: Path, sequence_min_length: int, prefix=None) -> Path:
    """
    Removes duplicates in no_phages list, writes all records in list to output.

    Verifies length of each sequence is longer than SEQUENCE_MIN_LENGTH.

    :param records: list of records from all protein files without keyword "phage"
    :param output: Path to output directory for profile HMMs and intermediate files
    :param prefix: Prefix for intermediate and result files
    :param sequence_min_length: Minimum sequence length for a record to be included in the input
    :return: Path to curated FASTA file without repeats or phages
    """
    output_dir = output / Path("intermediate_files")

    if not output_dir.exists():
        output_dir.mkdir()

    output_name = "curated_records.faa"
    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = output_dir / Path(output_name)
    SeqIO.write(remove_dupes(records, sequence_min_length), Path(output_path), "fasta")

    return output_path
