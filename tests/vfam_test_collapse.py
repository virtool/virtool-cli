import filecmp
import pytest
import os

from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_curate import *
from pathlib import Path

OUTPUT_PATH = "vfam_output"
TEST_PATH = "vfam_input"
CLUSTER_FILE_PATH = "vfam_intermediates/generic_input.clstr"
BLAST_FILE_PATH = "vfam_intermediates/generic_input_all_by_all_blast.br"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


def test_cluster_sequences():
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""
    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, Path(OUTPUT_PATH))

    cdhit_result = generate_clusters(no_dupes, Path(OUTPUT_PATH))

    assert filecmp.cmp(cdhit_result, Path(CLUSTER_FILE_PATH), shallow=True)






def test_all_by_all_blast():
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""
    pass


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v","-s"])