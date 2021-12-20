from collections import defaultdict
from pathlib import Path
from typing import List, Optional
from Bio import SearchIO
from virtool_cli.vfam_console import console


def get_sequence_lengths(blast_results_path: Path) -> dict:
    """
    Takes path to BLAST results file in BLAST tabular output format 6, parses file using SearchIO.

    If query sequence ID matches the subject sequence ID, sequence length is given as the length of the query_result.

    Sequence lengths of each sequence are stored as {sequence ID: sequence length} pairs in seq_lengths dictionary.

    :param blast_results_path: path to BLAST file produced in blast_all_by_all step
    :return: sequence_lengths, a dictionary containing each sequence ID and its sequence length
    """
    seq_lengths = dict()

    for query_result in SearchIO.parse(blast_results_path, "blast-tab"):
        for hit in query_result.hits:
            if query_result.id == hit.id:
                for hsp in hit.hsps:
                    seq_lengths[query_result.id] = hsp.aln_span

    return seq_lengths


def get_alignment_records(blast_results_path: Path) -> defaultdict:
    """
    Takes path to BLAST file in BLAST tabular output format 6 and parses file using SearchIO.

    Iterates through QueryResult objects to gather query ID, subject ID, query start, and query end info.

    QueryResult information for each query_result is stored in alignment_records dictionary

    :param blast_results_path: path to BLAST file produced in blast_all_by_all step
    :return: alignment_records, a dictionary containing all query_result objects for each query
    """
    alignment_records = defaultdict(list)

    for query_result in SearchIO.parse(blast_results_path, "blast-tab"):
        for hit in query_result.hits:
            if query_result.id != hit.id:
                for hsp in hit.hsps:
                    alignment_records[query_result.id].append(
                        {
                            "q_id": query_result.id,
                            "s_id": hit.id,
                            "q_start": hsp.query_start,
                            "q_end": hsp.query_end,
                        }
                    )

    return alignment_records


def check_alignments_by_length(
    seq_id: str, alignment_records: dict, seq_lengths: dict
) -> List[defaultdict]:
    """
    Iterates through alignment records for seq_id, adds alignment record to checked_alignments if:

     - The length of the subject is less than 70% the length of the query
     - The alignment spans more than 70% of the subject length

    :param seq_id: sequence ID for which to filter alignment records
    :param alignment_records: dictionary containing all alignment information for each seq_id
    :param seq_lengths: dictionary containing {sequence ID: sequence} length pairs
    :return: checked alignments, a list of alignment records to investigate further in check_alignments_by_position()
    """
    checked_alignments = list()

    for alignment in alignment_records[seq_id]:
        if seq_lengths[alignment["s_id"]] < 0.7 * seq_lengths[alignment["q_id"]]:
            subject_coverage = (
                float(abs(alignment["q_start"] - alignment["q_end"]))
                / seq_lengths[alignment["s_id"]]
            )

            if subject_coverage >= 0.7:
                checked_alignments.append(alignment)

    return checked_alignments


def check_alignments_by_position(
    seq_id: str, checked_by_length: List[dict], seq_lengths: dict
) -> Optional[str]:
    """
    Iterates through alignment records for sequence ID to be further investigated from check_alignments_by_length().

    If alignment queries collectively span more than 80% of sequence length, sequence is identified as polyprotein-like.

    :param seq_id: sequence ID for sequence to be investigated
    :param checked_by_length: alignments to be further investigated from check_alignment_records()
    :param seq_lengths: dictionary containing {sequence ID: sequence length} pairs
    :return: seq_id if found to be polyprotein-like
    """
    query_coverage = dict()

    for alignment in checked_by_length:
        for position in range(alignment["q_start"], alignment["q_end"]):
            query_coverage[position] = None

    if len(query_coverage) > 0.8 * (seq_lengths[seq_id]):
        return seq_id

    return None


def find_polyproteins(blast_results_path: Path) -> List[str]:
    """
    Sequences are filtered by sequence length and coverage to determine if they are polyprotein_ids/polyprotein-like.

    Sequences longer than 400 amino acids in length were identified as polyprotein or polyprotein-like if
    - at least 70% of the sequence length was covered by two or more other proteins in the sequence set.
    - these two or more other proteins were covered at least 80% by the longer sequence.

    :param blast_results_path: path to BLAST file produced in blast_all_by_all() step
    :return: polyprotein_ids, a list of sequence IDs to not include in output
    """
    seq_lengths = get_sequence_lengths(blast_results_path)
    alignment_records = get_alignment_records(blast_results_path)
    polyprotein_ids = []

    for seq_id in seq_lengths:
        if seq_lengths[seq_id] > 400 and seq_id in alignment_records:

            checked_by_length = check_alignments_by_length(
                seq_id, alignment_records, seq_lengths
            )

            if checked_by_length:
                checked_by_position = check_alignments_by_position(
                    seq_id, checked_by_length, seq_lengths
                )

                if checked_by_position:
                    polyprotein_ids.append(checked_by_position)

    if len(polyprotein_ids) > 0:
        console.print(
            f"✔ Filtered out {len(polyprotein_ids)} polyprotein-like records based on coverage.",
            style="green",
        )
    return polyprotein_ids
