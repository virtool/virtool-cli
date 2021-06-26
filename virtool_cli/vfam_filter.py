from Bio import SeqIO
from pathlib import Path

from typing import Optional, List

COVERAGE_HEUR_DICT = {
    0: 0.6,
    1: 0.65,
    2: 0.7,
    3: 0.75,
    4: 0.8,
    5: 0.85
}


def filter_file_on_coverage(fasta_path: Path) -> Optional[Path]:
    """
    Takes a FASTA file and a dictionary containing coverage heuristics information.

    Any sequences that don't subscribe to these heuristics are removed.

    :param fasta_path: path to FASTA file containing unfiltered protein sequence records
    :return: path to filtered FASTA file
    """
    record_lengths = {}
    lengths = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        record_lengths[record.id] = len(record.seq)
        lengths.append(len(record.seq))

    lengths.sort()

    upper = lengths[int(len(lengths) / 2)]
    lower = lengths[int((len(lengths) / 2) - 1)]
    median = float(upper + lower) / 2

    if len(lengths) % 2 != 0:
        median = float(lengths[int(len(lengths) / 2)])

    coverage_key = int(median / 100)

    coverage_threshold = (1.0 + max(COVERAGE_HEUR_DICT.values())) / 2

    if coverage_key in COVERAGE_HEUR_DICT:
        coverage_threshold = (1.0 + COVERAGE_HEUR_DICT[coverage_key]) / 2

    to_remove = [seq_id for seq_id, length in record_lengths.items()
                 if length < coverage_threshold * median or length * coverage_threshold > median]

    for seq_id in to_remove:
        if seq_id in record_lengths:
            record_lengths.pop(seq_id)

    if len(record_lengths) == 0:
        return None

    if len(record_lengths) == len(lengths):
        return fasta_path

    output_path = Path(f"{fasta_path}_filtered")
    to_write = (record for record in SeqIO.parse(fasta_path, "fasta") if record.id in record_lengths)
    SeqIO.write(to_write, output_path, "fasta")

    return output_path


def filter_on_coverage(fasta_paths: List[Path]) -> List[Path]:
    """
    Takes in list of FASTA paths, and calls filter_file_on_coverage to filter each file by coverage.

    :param fasta_paths: list of paths to FASTA files from mcl_to_fasta step to be filtered
    :return: filtered_by_coverage, a list of paths to filtered FASTA files
    """
    filtered_by_coverage = []
    for fasta_path in fasta_paths:
        filtered_file = filter_file_on_coverage(fasta_path)
        if filtered_file:
            filtered_by_coverage.append(filtered_file)

    return filtered_by_coverage


def filter_on_number(fasta_paths: List[Path], min_sequences: int) -> List[Path]:
    """
    Takes in a list of FASTA paths and filters out files if they contain less than MIN_SEQUENCES records.

    :param fasta_paths: list of FASTA files to be filtered
    :param min_sequences: Filter out clusters with fewer records than min_sequences_check
    :return: filtered_fasta_paths, a list of filtered FASTA files
    """
    filtered_fasta_paths = []
    for fasta_path in fasta_paths:

        for index, _ in enumerate(SeqIO.parse(fasta_path, "fasta")):
            if index + 1 >= min_sequences:
                filtered_fasta_paths.append(fasta_path)
                break

    return filtered_fasta_paths
