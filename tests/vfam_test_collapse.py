import filecmp
import pytest


from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_curate import *
from pathlib import Path


DUPES_INPUT = "vfam_input/Dupes"
GENERIC_INPUT = "vfam_input/Generic"
LARGE_INPUT = "vfam_input/Large"

COLLAPSED_DUPES = "vfam_og_intermediates/Dupes/collapsed_fasta"
COLLAPSED_GENERIC = "vfam_og_intermediates/Generic/collapsed_fasta"
COLLAPSED_LARGE = "vfam_og_intermediates/Large/collapsed_fasta"

BLAST_DUPES = "vfam_og_intermediates/Dupes/blast.br"
BLAST_GENERIC = "vfam_og_intermediates/Generic/blast.br"
BLAST_LARGE = "vfam_og_intermediates/Large/blast.br"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.mark.parametrize("input_dir, clustered_file", [(DUPES_INPUT, COLLAPSED_DUPES), (GENERIC_INPUT, COLLAPSED_GENERIC), (LARGE_INPUT, COLLAPSED_LARGE)])
def test_cluster_sequences(input_dir, clustered_file, output):
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, 1)

    result = generate_clusters(no_dupes, None, 1.0)

    assert filecmp.cmp(result, Path(clustered_file))

    result = str(result) + ".clstr"
    clustered_file = str(clustered_file) + ".clstr"
    assert filecmp.cmp(result, clustered_file)


@pytest.mark.parametrize("input_dir, blast_file", [(DUPES_INPUT, BLAST_DUPES), (GENERIC_INPUT, BLAST_GENERIC), (LARGE_INPUT, BLAST_LARGE)])
def test_all_by_all_blast(input_dir, blast_file, output):
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""

    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, 1)

    clustered_file = generate_clusters(no_dupes, None, 1.0)
    result = all_by_all_blast(clustered_file, 8)
    assert filecmp.cmp(result, Path(blast_file))


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])
