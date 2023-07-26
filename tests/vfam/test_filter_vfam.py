import os
import pytest
import filecmp
from Bio import SeqIO
from pathlib import Path

from fixtures_vfam import *
from paths_vfam import VFAM_INTERMEDIATES_PATH, DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT

DUPES_COVERAGE = VFAM_INTERMEDIATES_PATH / "Dupes" / "filtered_on_coverage"
GENERIC_COVERAGE = VFAM_INTERMEDIATES_PATH / "Generic" / "filtered_on_coverage"
LARGE_COVERAGE = VFAM_INTERMEDIATES_PATH / "Large" / "filtered_on_coverage"

DUPES_REMOVED_CVG = VFAM_INTERMEDIATES_PATH / "Dupes" / "removed_on_coverage"
GENERIC_REMOVED_CVG = VFAM_INTERMEDIATES_PATH / "Generic" / "removed_on_coverage"
LARGE_REMOVED_CVG = VFAM_INTERMEDIATES_PATH / "Large" / "removed_on_coverage"

DUPES_FILTERED_BY_NUM = VFAM_INTERMEDIATES_PATH / "Dupes" / "filtered_by_number"
GENERIC_FILTERED_BY_NUM = VFAM_INTERMEDIATES_PATH / "Generic" / "filtered_by_number"


@pytest.mark.parametrize(
    "input_dir, filtered_files",
    [
        (DUPES_INPUT, DUPES_COVERAGE),
        (GENERIC_INPUT, GENERIC_COVERAGE),
        (LARGE_INPUT, LARGE_COVERAGE),
    ],
)
def test_on_cvg_files(input_dir, filtered_files, filtered_on_cvg):
    """Test that filter_on_coverage filters the same files as original vfam."""
    with Path(filtered_files).open("r") as handle:
        expected = [line.split()[0] for line in handle]
    expected.sort()
    result = filtered_on_cvg
    assert len(expected) == len(result)


@pytest.mark.parametrize(
    "input_dir, removed_recs",
    [
        (DUPES_INPUT, DUPES_REMOVED_CVG),
        (GENERIC_INPUT, GENERIC_REMOVED_CVG),
        (LARGE_INPUT, LARGE_REMOVED_CVG),
    ],
)
def test_on_cvg_recs(input_dir, removed_recs, filtered_on_cvg):
    """Test that filter_on_coverage removes the same records as original vfam."""
    with Path(removed_recs).open("r") as handle:
        removed_records = [line.split()[0] for line in handle]

    result = filtered_on_cvg
    for file in result:
        for record in SeqIO.parse(file.open("r"), "fasta"):
            assert record.id not in removed_records


@pytest.mark.parametrize(
    "input_dir, filtered_files",
    [(DUPES_INPUT, DUPES_FILTERED_BY_NUM), (GENERIC_INPUT, GENERIC_FILTERED_BY_NUM)],
)
def test_on_number(input_dir, filtered_files, filtered_on_num):
    """Test that filter_on_number filters the same files as original vfam."""
    result = filtered_on_num
    result.sort()

    filtered_names = os.listdir(filtered_files)
    expected = []
    for name in filtered_names:
        expected.append(filtered_files / Path(name))
    expected.sort()

    assert len(expected) == len(result)
    for r_file, e_file in zip(result, expected):
        assert filecmp.cmp(r_file, e_file)
