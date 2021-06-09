from virtool_cli.vfam_curate import get_input_paths, remove_phages, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, polyprotein_name_check, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import blast_to_mcl, mcl_to_fasta
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number
from pathlib import Path


def run(src_path, output, sequence_min_length, phage_name_check, fraction_coverage, fraction_id, num_cores,
        polyp_name_check, inflation_num, filter_on_cvg, min_sequences):
    """
    Dictates workflow for vfam pipeline
    """
    # steps in vfam_curate module
    input_paths = get_input_paths(Path(src_path))

    if phage_name_check:
        no_phages = remove_phages(input_paths)
    else:
        no_phages = group_input_paths(input_paths)

    no_dupes = remove_dupes(no_phages, output, sequence_min_length)

    # steps found in vfam_collapse module
    cd_hit_result = generate_clusters(no_dupes, fraction_coverage, fraction_id)

    if polyp_name_check:
        cd_hit_result = polyprotein_name_check(cd_hit_result)

    blast_results = all_by_all_blast(cd_hit_result, num_cores)

    # steps found in vfam_polyprotein module
    polyproteins = find_polyproteins(blast_results)

    # steps found in vfam_markov module
    mcl_results = blast_to_mcl(blast_results, polyproteins, inflation_num)

    fasta_files = mcl_to_fasta(mcl_results, cd_hit_result)

    # steps found in vfam_filter module
    if filter_on_cvg:
        fasta_files = filter_on_coverage(fasta_files)

    fasta_files = filter_on_number(fasta_files, min_sequences)






