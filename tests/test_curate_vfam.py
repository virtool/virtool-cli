import pytest
import filecmp
import os

from Bio import SeqIO
from pathlib import Path
from vfam_paths import VFAM_INPUT_PATH, VFAM_INTERMEDIATES_PATH
from virtool_cli.vfam_curate import group_input_paths, write_curated_recs

DUPES_INPUT = VFAM_INPUT_PATH / "Dupes"
GENERIC_INPUT = VFAM_INPUT_PATH / "Generic"
LARGE_INPUT = VFAM_INPUT_PATH / "Large"

FILTERED_DUPES = VFAM_INTERMEDIATES_PATH / "Dupes" / "filtered_fasta"
FILTERED_GENERIC = VFAM_INTERMEDIATES_PATH / "Generic" / "filtered_fasta"
FILTERED_LARGE = VFAM_INTERMEDIATES_PATH / "Large" / "filtered_fasta"


@pytest.mark.parametrize("input_dir", [DUPES_INPUT,
                                       GENERIC_INPUT,
                                       LARGE_INPUT])
def test_get_input_paths(input_dir, input_paths):
    """Test that correct file names are found after calling get_input_paths."""
    listdir = os.listdir(input_dir)
    assert len(listdir) == len(input_paths)

    for path in input_paths:
        assert os.path.split(path)[1] in listdir


@pytest.mark.parametrize("input_dir", [DUPES_INPUT,
                                       GENERIC_INPUT,
                                       LARGE_INPUT])
def test_remove_phages(input_paths, input_dir):
    """Test that "phage" is not found in record descriptions for any records in records list."""
    no_phages = group_input_paths(input_paths, True)
    for record in no_phages:
        assert "phage" not in record.description


@pytest.mark.parametrize("input_dir", [DUPES_INPUT,
                                       GENERIC_INPUT,
                                       LARGE_INPUT])
def test_curated_file(input_dir, group_records, output):
    """Assert all sequences are longer than min_length, no record descriptions contain 'phage'."""
    no_dupes = write_curated_recs(group_records, output, None, 1)
    with Path(no_dupes) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            assert len(record.seq) > 1


@pytest.mark.parametrize("input_dir, filtered_file", [(DUPES_INPUT, FILTERED_DUPES),
                                                      (GENERIC_INPUT, FILTERED_GENERIC),
                                                      (LARGE_INPUT, FILTERED_LARGE)])
def test_remove_dupes(input_dir, input_paths, filtered_file, output):
    """Test that remove_dupes step of program produces desired output to match filtered og vfam file."""
    records = group_input_paths(input_paths, False)
    result = write_curated_recs(records, output, None, 1)
    expected = Path(filtered_file)

    assert filecmp.cmp(result, expected, shallow=True)
