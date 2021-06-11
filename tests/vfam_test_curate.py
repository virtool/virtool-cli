import pytest
import filecmp
import os


from Bio import SeqIO
from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, remove_phages, group_input_paths, remove_dupes


DUPES_INPUT = "vfam_input/Dupes"
GENERIC_INPUT = "vfam_input/Generic"
LARGE_INPUT = "vfam_input/Large"

FILTERED_DUPES = "vfam_og_intermediates/Dupes/filtered_fasta"
FILTERED_GENERIC = "vfam_og_intermediates/Generic/filtered_fasta"
FILTERED_LARGE = "vfam_og_intermediates/Large/filtered_fasta"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def input_paths(input_dir):
    return get_input_paths(Path(input_dir))


@pytest.fixture()
def no_phages(input_paths):
    return remove_phages(input_paths)


@pytest.fixture()
def no_dupes(no_phages, output):
    return remove_dupes(no_phages, output, None, 1)


@pytest.mark.parametrize("input_dir", [DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT])
def test_get_input_paths(input_dir, input_paths):
    """
    Test that same number of files and correct file names are found in input_paths after calling get_input_paths
    """
    listdir = os.listdir(input_dir)
    assert len(listdir) == len(input_paths)

    for path in input_paths:
        assert os.path.split(path)[1] in listdir


@pytest.mark.parametrize("input_dir", [DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT])
def test_remove_phages(input_dir, no_phages):
    """
    Test that "phage" is not found in record descriptions for any records in records list
    """
    for record in no_phages:
        assert "phage" not in record.description


@pytest.mark.parametrize("input_dir", [DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT])
def test_curated_file(input_dir, output, no_dupes):
    """
    Assert all sequences are longer than min_length

    Assert that no record descriptions contain "phage" when remove_phages is called
    """
    with Path(no_dupes) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            assert "phage" not in record.description
            assert len(record.seq) > 1


@pytest.mark.parametrize("input_dir, filtered_file", [(DUPES_INPUT, FILTERED_DUPES), (GENERIC_INPUT, FILTERED_GENERIC),
                                                      (LARGE_INPUT, FILTERED_LARGE)])
def test_remove_dupes(input_dir, filtered_file, output):
    """
    Test that remove_dupes step of program produces desired output to match filtere og vfam file
    """
    input_paths = get_input_paths(Path(input_dir))
    records = group_input_paths(input_paths)
    result = remove_dupes(records, output, None, 1)
    expected = Path(filtered_file)

    assert filecmp.cmp(result, expected, shallow=True)

