from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, write_curated_recs
from virtool_cli.vfam_collapse import generate_clusters, write_rmv_polyproteins, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import blast_to_mcl, mcl_to_fasta
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number
from virtool_cli.vfam_msa import batch_muscle_call, batch_hmm_call, concatenate_hmms, organize_intermediates
from virtool_cli.vfam_annotation import clusters_to_json


def run(src_path, output, prefix, sequence_min_length, no_named_phages, fraction_coverage, fraction_id, num_cores,
        no_named_polyproteins, inflation_num, coverage_check, min_sequences):
    """Dictates workflow for vfam pipeline."""
    input_paths = get_input_paths(Path(src_path))

    records = group_input_paths(input_paths, no_named_phages)

    no_dupes = write_curated_recs(records, output, prefix, sequence_min_length)

    cd_hit_result = generate_clusters(no_dupes, prefix, fraction_coverage, fraction_id)

    if no_named_polyproteins:
        cd_hit_result = write_rmv_polyproteins(cd_hit_result, prefix)

    blast_results = all_by_all_blast(cd_hit_result, prefix, num_cores)

    polyproteins = find_polyproteins(blast_results)

    mcl_results = blast_to_mcl(blast_results, polyproteins, inflation_num, prefix)

    fasta_files = mcl_to_fasta(mcl_results, cd_hit_result, prefix)

    if coverage_check:
        fasta_files = filter_on_coverage(fasta_files)

    fasta_files = filter_on_number(fasta_files, min_sequences)

    aligned_files = batch_muscle_call(fasta_files)

    hmm_files = batch_hmm_call(aligned_files)

    concatenate_hmms(hmm_files, output, prefix)

    clusters_to_json(fasta_files, output)

    organize_intermediates(output)







