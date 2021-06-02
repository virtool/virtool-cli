from virtool_cli.vfam_curate import *
from virtool_cli.vfam_collapse import *


def run(src_path, output, fraction_coverage, fraction_id):
    """
    Dictates workflow for vfam pipeline
    """
    # steps in vfam_curate module
    input_paths = get_input_paths(src_path)
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    # steps found in vfam_collapse module
    clustered_fasta = generate_clusters(no_dupes, output, fraction_coverage, fraction_id)
    blast_db = all_by_all_blast(clustered_fasta)


