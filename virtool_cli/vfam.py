from pathlib import Path

from virtool_cli.vfam_curate import get_input_paths, group_input_paths, write_curated_recs
from virtool_cli.vfam_collapse import generate_clusters, write_rmv_polyproteins, blast_all_by_all
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import blast_to_mcl, mcl_to_fasta
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number
from virtool_cli.vfam_msa import batch_muscle_call, batch_hmm_call, concatenate_hmms
from virtool_cli.vfam_annotation import get_json_from_clusters


def run(src_path, output, prefix, sequence_min_length, no_named_phages, fraction_coverage, fraction_id, num_cores,
        no_named_polyproteins, inflation_num, filter_clusters, min_sequences):
    """Dictates workflow for vfam pipeline."""
    input_paths = get_input_paths(Path(src_path))

    records = group_input_paths(input_paths, no_named_phages)

    no_dupes_path = write_curated_recs(records, output, sequence_min_length, prefix)

    cd_hit_result_path = generate_clusters(no_dupes_path, fraction_id, prefix, fraction_coverage)

    if no_named_polyproteins:
        cd_hit_result_path = write_rmv_polyproteins(cd_hit_result_path, prefix)

    blast_results_path = blast_all_by_all(cd_hit_result_path, num_cores, prefix)

    polyprotein_ids = find_polyproteins(blast_results_path)

    mcl_results_path = blast_to_mcl(blast_results_path, polyprotein_ids, inflation_num, prefix)

    fasta_paths = mcl_to_fasta(mcl_results_path, cd_hit_result_path, prefix)

    if filter_clusters:
        fasta_paths = filter_on_coverage(fasta_paths)

    fasta_paths = filter_on_number(fasta_paths, min_sequences)

    aligned_file_paths = batch_muscle_call(fasta_paths)

    hmm_file_paths = batch_hmm_call(aligned_file_paths)

    concatenate_hmms(hmm_file_paths, output, prefix)

    get_json_from_clusters(fasta_paths, output)







