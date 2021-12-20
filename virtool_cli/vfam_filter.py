from Bio import SeqIO
from pathlib import Path

from typing import Optional, List
from virtool_cli.vfam_console import console

COVERAGE_HEUR_DICT = {0: 0.6, 1: 0.65, 2: 0.7, 3: 0.75, 4: 0.8, 5: 0.85}


def remove_on_coverage(
    record_lengths: dict, median: float, coverage_threshold: float
) -> List[str]:
    """
    Filters records based on coverage threshold, returns record IDs that subscribe to coverage heuristics.

    :param record_lengths: dictionary of record ID, sequence length key-value pairs
    :param median: median sequence length from get_median()
    :param coverage_threshold: coverage threshold calculated in get_coverage_threshold()
    :return: list of record IDs for filtered records

    """
    filtered_record_ids = list()

    for seq_id, length in record_lengths.items():
        if (
            length >= coverage_threshold * median
            and length * coverage_threshold <= median
        ):
            filtered_record_ids.append(seq_id)

    return filtered_record_ids


def get_median(lengths: List[int]) -> float:
    """
    Calculates median sequence length from sequence lengths in lengths.

    :param lengths: list of sequence lengths
    :return: median
    """
    lengths.sort()

    upper = lengths[int(len(lengths) / 2)]
    lower = lengths[int((len(lengths) / 2) - 1)]
    median = float(upper + lower) / 2

    if len(lengths) % 2 != 0:
        median = float(lengths[int(len(lengths) / 2)])

    return median


def get_coverage_threshold(median: float) -> float:
    """
    Calculates coverage threshold from median and coverage heuristics dictionary.

    :param median:median sequence length calculated in get_median()
    :return: coverage threshold
    """
    coverage_key = int(median / 100)

    coverage_threshold = (1.0 + max(COVERAGE_HEUR_DICT.values())) / 2

    if coverage_key in COVERAGE_HEUR_DICT:
        coverage_threshold = (1.0 + COVERAGE_HEUR_DICT[coverage_key]) / 2

    return coverage_threshold


def filter_file_on_coverage(fasta_path: Path) -> Optional[Path]:
    """
    Takes a path to a FASTA file, records from file are filtered by median sequence length and coverage threshold.

    Filtered records are written to output_path.

    :param fasta_path: FASTA file containing unfiltered protein sequence records
    :return: output_path, the path to the filtered FASTA file
    """
    record_lengths = {}

    for record in SeqIO.parse(fasta_path, "fasta"):
        record_lengths[record.id] = len(record.seq)

    lengths = list(record_lengths.values())
    median = get_median(lengths)

    coverage_threshold = get_coverage_threshold(median)

    filtered_record_ids = remove_on_coverage(record_lengths, median, coverage_threshold)

    if filtered_record_ids:
        output_path = Path(f"{fasta_path}_filtered")
        to_write = (
            record
            for record in SeqIO.parse(fasta_path, "fasta")
            if record.id in filtered_record_ids
        )
        SeqIO.write(to_write, output_path, "fasta")

        return output_path

    return None


def filter_on_coverage(fasta_paths: List[Path]) -> List[Path]:
    """
    Takes in list of FASTA paths, and calls filter_file_on_coverage to filter each file by coverage.

    :param fasta_paths: list of paths to FASTA files from mcl_to_fasta step to be filtered
    :return: filtered_by_coverage, a list of paths to filtered FASTA files
    """
    num_unfiltered = len(fasta_paths)
    filtered_by_coverage = []
    for fasta_path in fasta_paths:
        filtered_file = filter_file_on_coverage(fasta_path)
        if filtered_file:
            filtered_by_coverage.append(filtered_file)

    num_filtered = len(filtered_by_coverage)

    console.print(
        f"✔ Filtered out {num_unfiltered - num_filtered} FASTA cluster files based on coverage.",
        style="green",
    )

    return filtered_by_coverage


def filter_on_number(fasta_paths: List[Path], min_sequences: int) -> List[Path]:
    """
    Takes in a list of FASTA paths and filters out files if they contain less than MIN_SEQUENCES records.

    :param fasta_paths: list of FASTA files to be filtered
    :param min_sequences: Filter out clusters with fewer records than min_sequences_check
    :return: filtered_fasta_paths, a list of filtered FASTA files
    """
    num_unfiltered = len(fasta_paths)
    filtered_fasta_paths = []
    for fasta_path in fasta_paths:

        for index, _ in enumerate(SeqIO.parse(fasta_path, "fasta")):
            if index + 1 >= min_sequences:
                filtered_fasta_paths.append(fasta_path)
                break

    num_filtered = len(filtered_fasta_paths)

    console.print(
        f"✔ Filtered out {num_unfiltered - num_filtered} "
        f"FASTA cluster files based on number of sequences.",
        style="green",
    )

    return filtered_fasta_paths
