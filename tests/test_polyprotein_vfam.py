import pytest


from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins


DUPES_INPUT = Path(__file__).parent.parent / "tests" / "vfam_input" / "Dupes"
GENERIC_INPUT = Path(__file__).parent.parent / "tests" / "vfam_input" / "Generic"
LARGE_INPUT = Path(__file__).parent.parent / "tests" / "vfam_input" / "Large"

POLYP_DUPES = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Dupes" / "polyproteins"
POLYP_GENERIC = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Generic" / "polyproteins"
POLYP_LARGE = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Large" / "polyproteins"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def polyproteins(input_dir, output):
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, None, 1)

    clustered_file = generate_clusters(no_dupes, None, None, 1.0)
    blast_file = all_by_all_blast(clustered_file, None, 8)
    return find_polyproteins(blast_file)


@pytest.mark.parametrize("input_dir, expected_polyproteins", [(DUPES_INPUT, POLYP_DUPES),
                                                              (GENERIC_INPUT, POLYP_GENERIC),
                                                              (LARGE_INPUT, POLYP_LARGE)])
def test_polyprotein_list(input_dir, expected_polyproteins, polyproteins, output):
    """Test that find_polyproteins catches same polyprotein sequences as original vfam."""
    result = polyproteins

    with Path(expected_polyproteins).open('r') as handle:
        expected = [line.split()[0] for line in handle]

    assert len(expected) == len(result)
    assert sorted(expected) == sorted(result)
