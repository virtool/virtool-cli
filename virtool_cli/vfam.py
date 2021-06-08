from virtool_cli.vfam_curate import *
from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_polyprotein import *
from virtool_cli.vfam_markov import *


def run(src_path, output, sequence_min_length, phage_name_check, fraction_coverage, fraction_id, num_cores, polyp_name_check, inflation_num):
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

    cluster_files = mcl_to_fasta(mcl_results, cd_hit_result)


