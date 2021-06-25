from pathlib import Path

from typing import List, Optional
from virtool_cli.vfam_polyprotein import Alignment


def get_sequence_lengths(blast_results: Path) -> dict:
    """
    Takes blast results file and parses through lines, creating an alignment object from each line

    If alignment query matches the alignment subject, sequence length is store print("Done!")d as the length of that query

    :param blast_results: blast file produced in all_by_all blast step
    :return: sequence_lengths, a dictionary containing each sequence and its sequence length
    """
    seq_lengths = {}

    with blast_results.open("r") as handle:
        for line in handle:
            alignment = Alignment(line)
            if alignment.query == alignment.subject:
                seq_lengths[alignment.query] = alignment.length

    return seq_lengths


def get_alignment_records(blast_results: Path) -> dict:
    """
    Takes blast file and parses through lines, producing an alignment object from each line

    If alignment query does not match subject, alignment is added to list of alignments for each query

    :param blast_results: blast file produced in all_by_all blast step
    :return: alignment_records, a dictionary containing all alignment objects for each query
    """
    alignment_records = {}
    with blast_results.open("r") as handle:

        for line in handle:
            alignment = Alignment(line)

            if alignment.query != alignment.subject:
                if alignment.query not in alignment_records:
                    alignment_records[alignment.query] = [alignment]
                else:
                    alignment_records[alignment.query].append(alignment)

    return alignment_records


def check_alignments_by_length(seq_id: str, alignment_records: dict, seq_lengths: dict) -> List[Alignment]:
    """
    Iterates through alignment records for seq_id, adds alignment record to checked_alignments if:

     - The length of the subject is less than 70% the length of the query
     - Alignment spans more than 70% of the alignment subject

    :param seq_id: sequence ID for which to filter alignment records
    :param alignment_records: dictionary containing all alignment records for each seq_id
    :param seq_lengths: dictionary containing sequence ID: sequence length pairs
    :return: checked alignments, a list of alignment records to investigate further in check_alignments_by_position()
    """
    checked_alignments = []

    for alignment in alignment_records[seq_id]:
        if seq_lengths[alignment.subject] < 0.7 * seq_lengths[alignment.query]:
            subject_coverage = float(abs(alignment.qstart - alignment.qend)) / seq_lengths[alignment.subject]

            if subject_coverage >= 0.7:
                checked_alignments.append(alignment)

    return checked_alignments


def check_alignments_by_position(seq_id: str, checked_by_length: List[Alignment], seq_lengths: dict) -> Optional[str]:
    """
    Iterates through alignment records for sequence ID to be further investigated from check_alignments_by_length().

    If alignment queries collectively span more than 80% of sequence length, sequence is identified as polyprotein-like.

    :param seq_id: sequence ID for sequence to be investigated
    :param checked_by_length: alignments to be further investigated from check_alignment_records()
    :param seq_lengths: dictionary containing sequence ID: sequence length pairs
    :return: seq_id if found polyprotein-like
    """
    query_coverage = dict()

    for alignment in checked_by_length:
        for position in range(alignment.qstart, alignment.qend):
            query_coverage[position] = None

    if len(query_coverage) > 0.8 * (seq_lengths[seq_id]):
        return seq_id

    return None


def find_polyproteins(blast_results_path: Path) -> List[str]:
    """
    Sequences longer than 400 amino acids are filtered by length and coverage to determine if they are polyprotein-like.

    :param blast_results_path: path to BLAST file produced in all_by_all blast step
    :return: polyprotein_ids, a list of sequence IDs for polyprotein-like records
    """
    seq_lengths = get_sequence_lengths(blast_results_path)
    alignment_records = get_alignment_records(blast_results_path)
    polyprotein_ids = []

    for seq_id in seq_lengths:
        if seq_lengths[seq_id] > 400 and seq_id in alignment_records:

            checked_by_length = check_alignments_by_length(seq_id, alignment_records, seq_lengths)

            if checked_by_length:
                checked_by_position = check_alignments_by_position(seq_id, checked_by_length, seq_lengths)

                if checked_by_position:
                    polyprotein_ids.append(checked_by_position)

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
