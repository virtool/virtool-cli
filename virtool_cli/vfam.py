from virtool_cli.vfam_curate import *
from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_polyprotein import *

def run(src_path, output, fraction_coverage, fraction_id, num_cores,polyp_name_check):
    """
    Dictates workflow for vfam pipeline
    """
    # steps in vfam_curate module
    input_paths = get_input_paths(Path(src_path))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    # steps found in vfam_collapse module
    cd_hit_result = generate_clusters(no_dupes, output, fraction_coverage, fraction_id)

    if polyp_name_check:
        cd_hit_result = polyprotein_name_check(cd_hit_result, output)
    blast_results = all_by_all_blast(cd_hit_result, output, num_cores)

    # steps found in vfam_polyprotein module
    polyproteins = find_polyproteins(blast_results)


