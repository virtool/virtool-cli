import pytest

from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_curate import *
from virtool_cli.vfam_polyprotein import *
from pathlib import Path


LARGE_INPUT_PATH = "vfam_large_input"
POLYPROTEIN_PATH = "vfam_intermediates/polyprotein_list"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def polyprotein_file(tmp_path):
    polyprotein_path = tmp_path / "dir"
    polyprotein_path.mkdir()
    return polyprotein_path / "polyprotein_file"


@pytest.fixture()
def blast_results(output):
    input_paths = get_input_paths(Path(LARGE_INPUT_PATH))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    # doesn't do filter by name step
    cdhit_result = generate_clusters(no_dupes, output, None, 1.0)
    return all_by_all_blast(cdhit_result, output, 8)


def test_polyprotein_list(output, blast_results):
    expected_polyproteins = []
    f = open(Path(POLYPROTEIN_PATH))
    for line in f:
        expected_polyproteins.append(line.split()[0])

    result_list = find_polyproteins(blast_results)
    assert len(expected_polyproteins) == len(result_list)
    assert sorted(expected_polyproteins) == sorted(result_list)


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])
