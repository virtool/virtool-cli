import filecmp
import pytest


from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_curate import *
from pathlib import Path


OUTPUT_PATH = "vfam_output"
TEST_PATH = "vfam_input"
CLUSTER_FILE_PATH = "vfam_intermediates/generic_input.clstr"
CLUSTER_FASTA_PATH = "vfam_intermediates/clustered_file.faa"
BLAST_FILE_PATH = "vfam_intermediates/generic_input_all_by_all_blast.br"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


def test_cluster_sequences(output):
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""
    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    cdhit_result = generate_clusters(no_dupes, output, None, 1.0)

    assert filecmp.cmp(cdhit_result, Path(CLUSTER_FASTA_PATH))
    cluster_file = str(cdhit_result) + ".clstr"
    assert filecmp.cmp(Path(cluster_file), Path(CLUSTER_FILE_PATH), shallow=True)


def test_all_by_all_blast(output):
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""

    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    cdhit_result = generate_clusters(no_dupes, output, None, 1.0)
    blast_results_file = all_by_all_blast(cdhit_result, output, 8)

    assert filecmp.cmp(blast_results_file, Path(BLAST_FILE_PATH))


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])
