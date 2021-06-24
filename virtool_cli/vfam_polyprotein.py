from collections import defaultdict
from pathlib import Path
from typing import List


def sequence_lengths(blast_results_path: Path) -> dict:
    """
    Takes path to blast results file and parses through lines, creating an alignment object from each line.

    If alignment query matches the alignment subject, sequence length is stored as the length of that query.

    :param blast_results_path: path to blast file produced in all_by_all blast step
    :return: sequence_lengths, a dictionary containing each sequence and its sequence length
    """
    seq_lengths = {}

    with blast_results_path.open("r") as handle:
        for line in handle:
            alignment = Alignment(line)
            if alignment.query == alignment.subject:
                seq_lengths[alignment.query] = alignment.length

    return seq_lengths


def get_alignment_records(blast_results_path: Path) -> dict:
    """
    Takes path to blast file and parses through lines, producing an alignment object from each line.

    If alignment query does not match subject, alignment is added to list of alignments for each query.

    :param blast_results_path: path to blast file produced in all_by_all blast step
    :return: alignment_records, a dictionary containing all alignment objects for each query
    """
    alignment_records = {}
    with blast_results.open("r") as handle:

        for line in handle:
            alignment = Alignment(line)

            if alignment.query != alignment.subject:
                alignment_records[alignment.query].append(alignment)

    return alignment_records


def find_polyproteins(blast_results_path: Path) -> List[str]:
    """
    Sequences are filtered by sequence length and coverage to determine if they are polyprotein_ids/polyprotein-like.

    Sequences longer than 400 amino acids in length were identified as polyprotein or polyprotein-like if

    - at least 70% of the sequence length was covered by two or more other proteins in the sequence set.
    - these two or more other proteins were covered at least 80% by the longer sequence.

    :param blast_results_path: path to blast file produced in all_by_all blast step
    :return: polyprotein_ids, a list of sequences to not include in output
    """
    seq_lengths = sequence_lengths(blast_results_path)
    alignment_records = get_alignment_records(blast_results_path)
    polyprotein_ids = []

    for query in seq_lengths:
        if seq_lengths[query] > 400 and query in alignment_records:

            alignment_ranges = []

            for alignment in alignment_records[query]:
                if alignment.subject in seq_lengths:
                    if seq_lengths[alignment.subject] < 0.7 * seq_lengths[alignment.query]:
                        subject_cvg = float(abs(alignment.qstart - alignment.qend))/seq_lengths[alignment.subject]
                        if subject_cvg >= 0.7:
                            alignment_ranges.append((alignment.qstart, alignment.qend))

            query_cvg = {}
            for rng in alignment_ranges:
                for position in range(rng[0], rng[1]):
                    query_cvg[position] = None

            if len(query_cvg) > 0.8 * (seq_lengths[query]):
                polyprotein_ids.append(query)
            query_cvg.clear()

    return polyprotein_ids


class Alignment:
    """
    Class to facilitate parsing blast tabular output format 6

    Naming conventions and descriptions from https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    """
    def __init__(self, blast_data):
        """Assigns field names to data gathered from lines in tab-delimited format."""
        blast_data = blast_data.split("\t")
        self.query = blast_data[0]
        self.subject = blast_data[1]
        self.pident = float(blast_data[2])
        self.length = float(blast_data[3])
        self.mismatch = float(blast_data[4])
        self.gapopen = float(blast_data[5])
        self.qstart = int(blast_data[6])
        self.qend = int(blast_data[7])
        self.sstart = int(blast_data[8])
        self.send = int(blast_data[9])
        self.evalue = str(blast_data[10])
        self.bitscore = float(blast_data[11])
