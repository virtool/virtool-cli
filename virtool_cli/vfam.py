from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, remove_phages, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, polyprotein_name_check, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import blast_to_mcl, mcl_to_fasta
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number
from virtool_cli.vfam_msa import batch_muscle_call, batch_hmm_call, concatenate_hmms, organize_intermediates


def run(src_path, output, prefix, sequence_min_length, phage_name_check, fraction_coverage, fraction_id, num_cores,
        polyp_name_check, inflation_num, filter_on_cvg, min_sequences):
    """Dictates workflow for vfam pipeline."""
    input_paths = get_input_paths(Path(src_path))

    if phage_name_check:
        records = remove_phages(input_paths)
    else:
        records = group_input_paths(input_paths)

    no_dupes = remove_dupes(records, output, prefix, sequence_min_length)

    cd_hit_result = generate_clusters(no_dupes, prefix, fraction_coverage, fraction_id)

    if polyp_name_check:
        cd_hit_result = polyprotein_name_check(cd_hit_result, prefix)

    blast_results = all_by_all_blast(cd_hit_result, prefix, num_cores)

    polyproteins = find_polyproteins(blast_results)

    mcl_results = blast_to_mcl(blast_results, polyproteins, inflation_num, prefix)

    fasta_files = mcl_to_fasta(mcl_results, cd_hit_result, prefix)

    if filter_on_cvg:
        fasta_files = filter_on_coverage(fasta_files)

    fasta_files = filter_on_number(fasta_files, min_sequences)

    aligned_files = batch_muscle_call(fasta_files)

    hmm_files = batch_hmm_call(aligned_files)

    concatenate_hmms(hmm_files, output, prefix)

    organize_intermediates(output)







