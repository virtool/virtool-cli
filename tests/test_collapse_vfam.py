import filecmp
import pytest

from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, all_by_all_blast
from vfam_paths import VFAM_INPUT_PATH, VFAM_INTERMEDIATES_PATH

DUPES_INPUT = VFAM_INPUT_PATH / "Dupes"
GENERIC_INPUT = VFAM_INPUT_PATH / "Generic"
LARGE_INPUT = VFAM_INPUT_PATH / "Large"

COLLAPSED_DUPES = VFAM_INTERMEDIATES_PATH / "Dupes" / "collapsed_fasta"
COLLAPSED_GENERIC = VFAM_INTERMEDIATES_PATH / "Generic" / "collapsed_fasta"
COLLAPSED_LARGE = VFAM_INTERMEDIATES_PATH / "Large" / "collapsed_fasta"

BLAST_DUPES = VFAM_INTERMEDIATES_PATH / "Dupes" / "blast.br"
BLAST_GENERIC = VFAM_INTERMEDIATES_PATH / "Generic" / "blast.br"
BLAST_LARGE = VFAM_INTERMEDIATES_PATH / "Large" / "blast.br"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def clustered_result(input_dir, output):
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, None, 1)

    return generate_clusters(no_dupes, None, None, 1.0)


@pytest.fixture()
def blast_result(clustered_result):
    return all_by_all_blast(clustered_result, None, 8)


@pytest.mark.parametrize("input_dir, clustered_file", [(DUPES_INPUT, COLLAPSED_DUPES),
                                                       (GENERIC_INPUT, COLLAPSED_GENERIC),
                                                       (LARGE_INPUT, COLLAPSED_LARGE)])
def test_cluster_sequences(input_dir, clustered_file, clustered_result, output):
    """Test that cluster output file from cluster sequences matches desired output from original vfam."""
    result = clustered_result
    assert filecmp.cmp(result, Path(clustered_file))

    result = str(result) + ".clstr"
    clustered_file = str(clustered_file) + ".clstr"
    assert filecmp.cmp(result, clustered_file)


@pytest.mark.parametrize("input_dir, blast_file", [(DUPES_INPUT, BLAST_DUPES),
                                                   (GENERIC_INPUT, BLAST_GENERIC),
                                                   (LARGE_INPUT, BLAST_LARGE)])
def test_all_by_all_blast(input_dir, blast_file, blast_result, output):
    """Test that cluster output file from cluster sequences matches desired output from original vfam."""
    assert filecmp.cmp(blast_result, Path(blast_file))
