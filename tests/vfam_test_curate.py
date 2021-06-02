import pytest
import subprocess


from virtool_cli.vfam_curate import *
from pathlib import Path


TEST_PATH = "vfam_input"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


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


def test_curated_file(output):
    """
    Test that all sequences are longer than min_length and that no record descriptions contain "phage" in output
    """
    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    with Path(no_dupes) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            assert "phage" not in record.description
            assert len(record.seq) > 1










