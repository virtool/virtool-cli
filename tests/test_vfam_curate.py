import pytest
import subprocess


from virtool_cli.vfam_curate import *
from pathlib import Path


TEST_PATH = "../tests/vfam_files"
OUTPUT_TEST_PATH = "../tests/vfam_files/output"


def test_get_input_paths():
    """
    Test that same number of files and correct file names are found in input_paths after calling get_input_paths
    """
    listdir = os.listdir(TEST_PATH)
    input_paths = get_input_paths(Path(TEST_PATH))
    assert len(listdir) == len(input_paths)

    for path in input_paths:
        assert os.path.split(path)[1] in listdir


def test_remove_phages():
    """
    Test that "phage" is not found in record descriptions for any records in no_phages list
    """
    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)

    for record in no_phages:
        assert "phage" not in record.description


def test_curated_file():
    """
    Test that no files longer than min_length or no files with keyword "phage" are found in curated fasta file
    """

    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)
    remove_dupes(no_phages, Path(OUTPUT_TEST_PATH))

    with Path(OUTPUT_TEST_PATH) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            assert "phage" not in record.description
            assert len(record.seq) > 1











