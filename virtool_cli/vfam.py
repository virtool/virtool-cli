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
    # steps in vfam_curate module
    print("Gathering input")
    input_paths = get_input_paths(Path(src_path))

    if phage_name_check:
        print("Filtering phages by name")
        records = remove_phages(input_paths)
    else:
        records = group_input_paths(input_paths)

    print("Removing duplicate records")
    no_dupes = remove_dupes(records, output, prefix, sequence_min_length)

    # steps found in vfam_collapse module
    cd_hit_result = generate_clusters(no_dupes, prefix, fraction_coverage, fraction_id)

    if polyp_name_check:
        print("Filtering out polyproteins by name")
        cd_hit_result = polyprotein_name_check(cd_hit_result, prefix)

    blast_results = all_by_all_blast(cd_hit_result, prefix, num_cores)

    # steps found in vfam_polyprotein module
    print("Filtering out polyproteins by sequence length and sequence coverage")
    polyproteins = find_polyproteins(blast_results)

    # steps found in vfam_markov module
    print("Running MUSCLE on blast results")
    mcl_results = blast_to_mcl(blast_results, polyproteins, inflation_num, prefix)

    print("Splitting FASTA file based on MUSCLE results")
    fasta_files = mcl_to_fasta(mcl_results, cd_hit_result, prefix)

    # steps found in vfam_filter module
    if filter_on_cvg:
        fasta_files = filter_on_coverage(fasta_files)

    print("Filtering cluster FASTA files")
    fasta_files = filter_on_number(fasta_files, min_sequences)

    # steps found in vfam_msa module
    print("Running MUSCLE on all cluster FASTA files")
    aligned_files = batch_muscle_call(fasta_files)

    print("Running hmmbuild on all alignment files")
    hmm_files = batch_hmm_call(aligned_files)

    print("Building master hmm file")
    hmm_file = concatenate_hmms(hmm_files, output, prefix)

    print(f"{hmm_file} contains profile-HMMS from all input clusters")

    print("Organising intermediate files")
    organize_intermediates(output)

    print("Done!")







