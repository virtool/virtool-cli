import os
import pytest
import filecmp

from pathlib import Path
from paths import VFAM_INPUT_PATH, VFAM_INTERMEDIATES_PATH

DUPES_INPUT = VFAM_INPUT_PATH / "Dupes"
GENERIC_INPUT = VFAM_INPUT_PATH / "Generic"
LARGE_INPUT = VFAM_INPUT_PATH / "Large"

DUPES_MSA = VFAM_INTERMEDIATES_PATH / "Dupes" / "msa_files"
GENERIC_MSA = VFAM_INTERMEDIATES_PATH / "Generic" / "msa_files"

DUPES_HMM = VFAM_INTERMEDIATES_PATH / "Dupes" / "hmm_files"
GENERIC_HMM = VFAM_INTERMEDIATES_PATH / "Generic" / "hmm_files"


@pytest.mark.parametrize(
    "input_dir, msa", [(DUPES_INPUT, DUPES_MSA), (GENERIC_INPUT, GENERIC_MSA)]
)
def test_batch_muscle_call(get_msa, input_dir, msa):
    """Test msa files by comparing to og vfam msa data."""
    msa_names = os.listdir(msa)
    expected = []
    for name in msa_names:
        expected.append(msa / Path(name))
    expected.sort()

    result = get_msa
    result.sort()

    assert len(expected) == len(result)
    for r_file, e_file in zip(result, expected):
        assert filecmp.cmp(r_file, e_file)


@pytest.mark.parametrize(
    "input_dir, hmm", [(DUPES_INPUT, DUPES_HMM), (GENERIC_INPUT, GENERIC_HMM)]
)
def test_batch_hmm_call(get_hmm, input_dir, hmm):
    """
    Test profile HMMs by comparing lines 17 onward (part that contains hmm data) 
    of output data and og vfam data.
    """
    hmm_names = os.listdir(hmm)
    expected = []
    for name in hmm_names:
        expected.append(hmm / Path(name))
    expected.sort()

    result = get_hmm
    result.sort()

    assert len(expected) == len(result)
    for r_file, e_file in zip(result, expected):
        with open(r_file, "r") as r_handle:
            r_lines = r_handle.readlines()
            r_lines = r_lines[16:]
            with open(e_file, "r") as e_handle:
                e_lines = e_handle.readlines()
                e_lines = e_lines[16:]

                assert r_lines == e_lines
