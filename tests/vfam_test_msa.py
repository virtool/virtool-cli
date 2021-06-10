import os
import pytest
import filecmp


from Bio import SeqIO
from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import blast_to_mcl, mcl_to_fasta
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number
from virtool_cli.vfam_msa import batch_muscle_call, batch_hmm_call, concatenate_hmms


DUPES_INPUT = "vfam_input/Dupes"
GENERIC_INPUT = "vfam_input/Generic"
LARGE_INPUT = "vfam_input/Large"


DUPES_MSA = "vfam_og_intermediates/Dupes/msa_files"
GENERIC_MSA = "vfam_og_intermediates/Generic/msa_files"

DUPES_HMM = "vfam_og_intermediates/Dupes/hmm_files"
GENERIC_HMM = "vfam_og_intermediates/Generic/hmm_files"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def get_msa(input_dir, output):
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, None, 1)
    clustered_file = generate_clusters(no_dupes, None, None, 1.0)
    blast_results = all_by_all_blast(clustered_file, None, 8)
    polyproteins = find_polyproteins(blast_results)
    mcl_file = blast_to_mcl(blast_results, polyproteins, None, None)
    fasta_files = mcl_to_fasta(mcl_file, clustered_file, None)
    fasta_files = filter_on_coverage(fasta_files)
    fasta_files = filter_on_number(fasta_files, 2)

    return batch_muscle_call(fasta_files)


@pytest.mark.parametrize("input_dir, MSA", [(DUPES_INPUT, DUPES_MSA), (GENERIC_INPUT, GENERIC_MSA)])
def test_batch_muscle_call(get_msa, input_dir, MSA):
    msa_names = os.listdir(MSA)
    expected = []
    for name in msa_names:
        expected.append(MSA / Path(name))
    expected.sort()

    result = get_msa
    result.sort()

    assert len(expected) == len(result)
    for x in range(len(expected)):
        assert filecmp.cmp(result[x], expected[x], shallow=True)


@pytest.mark.parametrize("input_dir, HMM", [(DUPES_INPUT, DUPES_HMM), (GENERIC_INPUT, GENERIC_HMM)])
def test_batch_hmm_call(get_msa, input_dir, HMM):
    """Test profile HMMs by comparing lines 17 onward (part that contains HMM data) of output data and og vfam data"""
    hmm_names = os.listdir(HMM)
    expected = []
    for name in hmm_names:
        expected.append(HMM / Path(name))
    expected.sort()

    result = batch_hmm_call(get_msa)
    result.sort()

    assert len(expected) == len(result)
    for x in range(len(expected)):
        with open(result[x], "r") as r_file:
            r_lines = r_file.readlines()
            r_lines = r_lines[16:]
            with open(expected[x], "r") as e_file:
                e_lines = e_file.readlines()
                e_lines = e_lines[16:]

                assert r_lines == e_lines


def test_concatenate_hmms():
    pass

