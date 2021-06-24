from Bio import SeqIO
from pathlib import Path

from typing import Optional, List

COVERAGE_HEUR_DICT = {0: 0.6, 1: 0.65, 2: 0.7, 3: 0.75, 4: 0.8, 5: 0.85}


def filter_file_on_coverage(fasta_file: Path) -> Optional[Path]:
    """
    Takes a fasta file and a dictionary containing coverage heuristics information

    any sequences that don't subscribe to these heuristics are removed

    :param fasta_file: fasta file containing unfiltered protein sequence records
    :return: filtered_fasta_file
    """
    record_lengths = {}
    lengths = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        record_lengths[record.id] = len(record.seq)
        lengths.append(len(record.seq))

    lengths.sort()

    upper = lengths[int(len(lengths) / 2)]
    lower = lengths[int((len(lengths) / 2) - 1)]
    median = float(upper + lower) / 2

    if len(lengths) % 2 != 0:
        median = float(lengths[int(len(lengths) / 2)])

    coverage_key = int(median / 100)

    coverage_threshold = (1.0 + max(COVERAGE_HEUR_DICT.values()))/2

    if coverage_key in COVERAGE_HEUR_DICT:
        coverage_threshold = (1.0 + COVERAGE_HEUR_DICT[coverage_key])/2

    to_remove = [ID for ID, length in record_lengths.items()
                 if length < coverage_threshold * median or length * coverage_threshold > median]

    for ID in to_remove:
        if ID in record_lengths:
            record_lengths.pop(ID)

    if len(record_lengths) == 0:
        return None
    elif len(record_lengths) == len(lengths):
        return fasta_file
    else:
        output_path = Path(f"{fasta_file}_filtered")
        to_write = [record for record in SeqIO.parse(fasta_file, "fasta") if record.id in record_lengths]
        SeqIO.write(to_write, output_path, "fasta")
        return output_path


def filter_on_coverage(fasta_files: list) -> List[Path]:
    """
    Takes in list of fasta files, and calls filter_file_on_coverage to filter them by coverage

    :param fasta_files: list of fasta_files from mcl_to_fasta step to be filtered
    :return: filtered_by_coverage, a list filtered fasta files
    """
    filtered_by_coverage = []
    for fasta_file in fasta_files:
        filtered_file = filter_file_on_coverage(fasta_file)
        if filtered_file:
            filtered_by_coverage.append(filtered_file)

    return filtered_by_coverage


def filter_on_number(fasta_files: List[Path], min_sequences: int) -> List[Path]:
    """
    Takes in a list of fasta files and filters out files if they contain less than MIN_SEQUENCES records

    :param fasta_files: list of fasta files to be filtered
    :param min_sequences: Filter out clusters with fewer records than min_sequences_check
    :return: filtered_files, a list of filtered fasta files
    """
    filtered_files = []
    for fasta_file in fasta_files:
        count = 0
        for _ in SeqIO.parse(fasta_file, "fasta"):
            count += 1
            if count >= min_sequences:
                filtered_files.append(fasta_file)

                break

    return filtered_files
