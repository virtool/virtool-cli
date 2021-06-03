def get_sequence_lengths(all_by_all_blast_results_file) -> dict:
    """
    Takes blast results file and parses through lines

    If alignment record found where subject is the same as query, sequence length is given by alignment length

    :param all_by_all_blast_results_file: blast file produced in all_by_all blast step
    :return: sequence_lengths, a dictionary containing each sequence and its sequence length
    """
    blast_file = open(all_by_all_blast_results_file)
    sequence_lengths = {}

    for line in blast_file:
        blast_data = line.split("\t")

        query = blast_data[0]
        subject = blast_data[1]
        sequence_alignment_length = int(blast_data[3])

        if query == subject:
            sequence_lengths[query] = sequence_alignment_length

    return sequence_lengths


def get_alignment_records(all_by_all_blast_results_file) -> dict:
    """
    Takes blast file and parses through lines

    If alignment record found where subject does not equal query, alignment is added to query key in alignment_records

    :param all_by_all_blast_results_file: blast file produced in all_by_all blast step
    :return: alignment_records, a dictionary containing all alignment records for each query
    """
    blast_file = open(all_by_all_blast_results_file)
    alignment_records = {}

    for line in blast_file:
        blast_data = line.split("\t")

        query = blast_data[0]
        subject = blast_data[1]
        start_of_query_alignment = int(blast_data[6])
        end_of_query_alignment = int(blast_data[7])

        if query != subject:
            if query not in alignment_records:
                alignment_records[query] = [(subject, start_of_query_alignment, end_of_query_alignment)]
            else:
                alignment_records[query].append((subject, start_of_query_alignment, end_of_query_alignment))

    return alignment_records


def get_polyproteins(all_by_all_blast_results_file) -> list:
    """
    Sequences longer than 400 amino acids in length were identified as polyprotein or polyprotein-like if

    - at least 70% of the sequence length was covered by two or more other proteins in the sequence set
    - these two or more other proteins were covered at least 80% by the longer sequence.

    :param all_by_all_blast_results_file:blast file produced in all_by_all blast step
    :return: polyprotein_sequences, a list of sequences to not include in output
    """
    sequence_lengths = get_sequence_lengths(all_by_all_blast_results_file)
    alignment_records = get_alignment_records(all_by_all_blast_results_file)

    polyprotein_sequences = []
    for query in sequence_lengths:
        if int(sequence_lengths[query]) > 400 and query in alignment_records:

            if len(alignment_records[query]) > 1:
                alignment_ranges = []

                for alignment in alignment_records[query]:
                    subject = alignment[0]
                    start_of_query_alignment = alignment[1]
                    end_of_query_alignment = alignment[2]

                    if subject in sequence_lengths:
                        if float(sequence_lengths[subject]) < 0.7 * float(sequence_lengths[query]):
                            subject_coverage = float(abs(start_of_query_alignment - end_of_query_alignment)) / \
                                               float(sequence_lengths[subject])

                            if subject_coverage >= 0.7:
                                alignment_ranges.append((start_of_query_alignment, end_of_query_alignment))

                query_coverage = {}
                for query_range in alignment_ranges:
                    for amino_acid_position in range(query_range[0], query_range[1]):
                        query_coverage[amino_acid_position] = None

                if len(query_coverage) > 0.8 * float(sequence_lengths[query]):
                    polyprotein_sequences.append(query)

                query_coverage.clear()

    return polyprotein_sequences
