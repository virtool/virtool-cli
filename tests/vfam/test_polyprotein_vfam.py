import pytest
from pathlib import Path

from vfam_fixtures import *
from paths import VFAM_INTERMEDIATES_PATH, DUPES_INPUT, GENERIC_INPUT, LARGE_INPUT

POLYP_DUPES = VFAM_INTERMEDIATES_PATH / "Dupes" / "polyprotein_ids"
POLYP_GENERIC = VFAM_INTERMEDIATES_PATH / "Generic" / "polyprotein_ids"
POLYP_LARGE = VFAM_INTERMEDIATES_PATH / "Large" / "polyprotein_ids"


@pytest.mark.parametrize(
    "input_dir, expected_polyproteins",
    [
        (DUPES_INPUT, POLYP_DUPES),
        (GENERIC_INPUT, POLYP_GENERIC),
        (LARGE_INPUT, POLYP_LARGE),
    ],
)
def test_polyprotein_list(input_dir, expected_polyproteins, polyproteins, output):
    """Test that find_polyproteins catches same polyprotein sequences as original vfam."""
    result = polyproteins

    with Path(expected_polyproteins).open("r") as handle:
        expected = [line.split()[0] for line in handle]

    assert len(expected) == len(result)
    assert sorted(expected) == sorted(result)
