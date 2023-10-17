import filecmp
import pytest

from pathlib import Path
from paths import VFAM_INPUT_PATH, VFAM_INTERMEDIATES_PATH

DUPES_INPUT = VFAM_INPUT_PATH / "Dupes"
GENERIC_INPUT = VFAM_INPUT_PATH / "Generic"
LARGE_INPUT = VFAM_INPUT_PATH / "Large"

COLLAPSED_DUPES = VFAM_INTERMEDIATES_PATH / "Dupes" / "collapsed_fasta"
COLLAPSED_GENERIC = VFAM_INTERMEDIATES_PATH / "Generic" / "collapsed_fasta"
COLLAPSED_LARGE = VFAM_INTERMEDIATES_PATH / "Large" / "collapsed_fasta"

BLAST_DUPES = VFAM_INTERMEDIATES_PATH / "Dupes" / "blast.br"
BLAST_GENERIC = VFAM_INTERMEDIATES_PATH / "Generic" / "blast.br"
BLAST_LARGE = VFAM_INTERMEDIATES_PATH / "Large" / "blast.br"


@pytest.mark.parametrize(
    "input_dir, clustered_file",
    [
        (DUPES_INPUT, COLLAPSED_DUPES),
        (GENERIC_INPUT, COLLAPSED_GENERIC),
        (LARGE_INPUT, COLLAPSED_LARGE),
    ],
)
def test_cluster_sequences(input_dir, input_paths, clustered_file, clustered_result):
    """Test that cluster output file from cluster sequences matches output from original vfam. Skip no_phages."""
    result = clustered_result
    assert filecmp.cmp(result, Path(clustered_file))

    result = str(result) + ".clstr"
    clustered_file = str(clustered_file) + ".clstr"
    assert filecmp.cmp(result, clustered_file)


@pytest.mark.parametrize(
    "input_dir, blast_file",
    [
        (DUPES_INPUT, BLAST_DUPES),
        (GENERIC_INPUT, BLAST_GENERIC),
        (LARGE_INPUT, BLAST_LARGE),
    ],
)
def test_all_by_all_blast(input_dir, input_paths, blast_result, blast_file):
    """Test that cluster output file from cluster sequences matches output from original vfam. Skip no_phages."""
    assert filecmp.cmp(blast_result, Path(blast_file))
