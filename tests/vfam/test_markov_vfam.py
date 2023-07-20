import pytest
import filecmp
import os
from pathlib import Path

from vfam_fixtures import *
from paths_vfam import VFAM_INTERMEDIATES_PATH, DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT

BLAST_DUPES = VFAM_INTERMEDIATES_PATH / "Dupes" / "blast"
BLAST_GENERIC = VFAM_INTERMEDIATES_PATH / "Generic" / "blast"
BLAST_LARGE = VFAM_INTERMEDIATES_PATH / "Large" / "blast"

DUPES_CLUSTERS = VFAM_INTERMEDIATES_PATH / "Dupes" / "fasta_files"
GENERIC_CLUSTERS = VFAM_INTERMEDIATES_PATH / "Generic" / "fasta_files"
LARGE_CLUSTERS = VFAM_INTERMEDIATES_PATH / "Large" / "fasta_files"


@pytest.mark.parametrize(
    "input_dir, blast_path",
    [
        (DUPES_INPUT, BLAST_DUPES),
        (GENERIC_INPUT, BLAST_GENERIC),
        (LARGE_INPUT, BLAST_LARGE),
    ],
)
def test_abc(input_dir, blast_path, output, abc_file):
    """Test that write_abc produces an abc files to match original vfam output."""
    result = Path(str(blast_path) + ".abc")

    assert filecmp.cmp(result, abc_file)


@pytest.mark.parametrize(
    "input_dir, blast_path",
    [
        (DUPES_INPUT, BLAST_DUPES),
        (GENERIC_INPUT, BLAST_GENERIC),
        (LARGE_INPUT, BLAST_LARGE),
    ],
)
def test_blast_to_mcl(input_dir, blast_path, output, mcl_results):
    """Test that blast_to_mcl produces .mcl files that match original vfam output."""
    result = mcl_results
    expected = str(blast_path) + ".mcl"
    assert filecmp.cmp(result, expected, shallow=True)


@pytest.mark.parametrize(
    "input_dir, blast_path",
    [
        (DUPES_INPUT, BLAST_DUPES),
        (GENERIC_INPUT, BLAST_GENERIC),
        (LARGE_INPUT, BLAST_LARGE),
    ],
)
def test_mci(input_dir, blast_path, output, mcl_results):
    """Test that blast_to_mcl produces .mci files that match original vfam output."""
    result = mcl_results.parent / "blast.mci"
    expected = str(blast_path) + ".mci"
    assert filecmp.cmp(result, expected, shallow=True)


@pytest.mark.parametrize(
    "input_dir, cluster_dir",
    [(DUPES_INPUT, DUPES_CLUSTERS), (GENERIC_INPUT, GENERIC_CLUSTERS)],
)
def test_mcl_to_fasta(input_dir, cluster_dir, output, mcl_results, fasta_files):
    """Test that mcl_to_fasta produces fasta files that match original vfam output."""
    result_files = fasta_files
    result_files.sort()

    expected_files = []
    for file in os.listdir(cluster_dir):
        expected_files.append(cluster_dir / Path(file))
    expected_files.sort()

    assert len(expected_files) == len(result_files)
    for r_file, e_file in zip(result_files, expected_files):
        assert filecmp.cmp(r_file, e_file)
