import pytest
import filecmp

from virtool_cli.vfam_curate import *
from pathlib import Path


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


@pytest.mark.parametrize("input_dir", [DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT])
def test_get_input_paths(input_dir):
    """
    Test that same number of files and correct file names are found in input_paths after calling get_input_paths
    """
    listdir = os.listdir(input_dir)
    input_paths = get_input_paths(Path(input_dir))
    assert len(listdir) == len(input_paths)

    for path in input_paths:
        assert os.path.split(path)[1] in listdir


@pytest.mark.parametrize("input_dir", [DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT])
def test_remove_phages(input_dir):
    """
    Test that "phage" is not found in record descriptions for any records in no_phages list
    """
    input_paths = get_input_paths(Path(input_dir))
    no_phages = remove_phages(input_paths)

    for record in no_phages:
        assert "phage" not in record.description


@pytest.mark.parametrize("input_dir", [DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT])
def test_curated_file(input_dir, output):
    """
    Test that all sequences are longer than min_length and that no record descriptions contain "phage" in output
    """
    input_paths = get_input_paths(Path(input_dir))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output, 1)

    with Path(no_dupes) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            assert "phage" not in record.description
            assert len(record.seq) > 1


@pytest.mark.parametrize("input_dir, filtered_file", [(DUPES_INPUT, FILTERED_DUPES), (GENERIC_INPUT, FILTERED_GENERIC), (LARGE_INPUT, FILTERED_LARGE)])
def test_remove_dupes(input_dir, filtered_file, output):
    """
    TODO: figure out why filecmp fails for dupes file
    Test that remove_dupes step of program produces desired output
    """
    input_paths = get_input_paths(Path(input_dir))
    records = group_input_paths(input_paths)
    result = remove_dupes(records, output, 1)
    expected = Path(filtered_file)

    result_phage_count = 0
    result_record_count = 0

    with result as handle:
        for record in SeqIO.parse(handle, "fasta"):
            result_record_count += 1
            if "phage" in record.description:
                result_phage_count += 1

    expected_phage_count = 0
    expected_record_count = 0
    with expected as handle:
        for record in SeqIO.parse(handle, "fasta"):
            expected_record_count += 1
            if "phage" in record.description:
                expected_phage_count += 1

    assert result_phage_count == expected_phage_count
    assert result_record_count == expected_record_count
    assert filecmp.cmp(result, expected)


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])
