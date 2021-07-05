from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, write_curated_recs, remove_dupes, \
    get_taxonomy, write_no_dupes
from virtool_cli.vfam_collapse import generate_clusters, write_rmv_polyproteins, blast_all_by_all
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import blast_to_mcl, mcl_to_fasta
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number
from virtool_cli.vfam_msa import batch_muscle_call, batch_hmm_call, concatenate_hmms
from virtool_cli.vfam_annotation import get_json_from_clusters


def run(src_path,
        output,
        prefix,
        sequence_min_length,
        no_named_phages,
        fraction_coverage,
        fraction_id,
        num_cores,
        no_named_polyproteins,
        inflation_num,
        coverage_check,
        min_sequences):
    """Dictates workflow for vfam pipeline."""
    input_paths = get_input_paths(Path(src_path))

    records = group_input_paths(input_paths, no_named_phages)

    no_dupes = remove_dupes(records, sequence_min_length)

    no_dupes_path = write_no_dupes(no_dupes, output, prefix)

    taxonomy_records = get_taxonomy(no_dupes_path)

    curated_recs = write_curated_recs(no_dupes_path, list(taxonomy_records.keys()), prefix)

    cd_hit_result = generate_clusters(curated_recs, fraction_id, fraction_coverage, prefix)

    if no_named_polyproteins:
        cd_hit_result = write_rmv_polyproteins(cd_hit_result, prefix)

    blast_results = blast_all_by_all(cd_hit_result, num_cores, prefix)

    polyproteins = find_polyproteins(blast_results)

    mcl_results = blast_to_mcl(blast_results, polyproteins, inflation_num, prefix)

    fasta_files = mcl_to_fasta(mcl_results, cd_hit_result, prefix)

    if coverage_check:
        fasta_files = filter_on_coverage(fasta_files)

    fasta_files = filter_on_number(fasta_files, min_sequences)

    aligned_files = batch_muscle_call(fasta_files)

    hmm_files = batch_hmm_call(aligned_files)

    concatenate_hmms(hmm_files, output, prefix)

    get_json_from_clusters(fasta_files, taxonomy_records, output)
