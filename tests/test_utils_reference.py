import pytest
from pathlib import Path

from virtool_cli.utils.reference import (
    is_v1,
    search_otu_by_id,
    get_otu_paths,
    get_isolate_paths,
)
from virtool_cli.utils.id_generator import generate_unique_ids
from paths import TEST_FILES_PATH

SRC_PATH = TEST_FILES_PATH / "src_test"
SRC_V1_PATH = TEST_FILES_PATH / "src_v1"


def get_subdir(path: Path) -> list:
    """
    Generates a list of all subdirectories under a path

    :param path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    return [subdir for subdir in path.iterdir() if subdir.is_dir()]


def test_utils_is_v1():
    assert is_v1(SRC_V1_PATH)


def test_utils_is_v2():
    assert not is_v1(SRC_PATH)


@pytest.mark.parametrize("otu_id", ["f8a56910", "3962f6ec"])
def test_utils_search_by_id_success(otu_id):
    assert search_otu_by_id(otu_id, SRC_PATH) is not None


@pytest.mark.parametrize("otu_id", ["bbbbbbbb", "777777tt"])
def test_utils_search_by_id_fail(otu_id):
    assert search_otu_by_id(otu_id, SRC_PATH) is None


@pytest.mark.parametrize("src", [SRC_PATH])
def test_utils_get_otu_paths(src):
    assert set(get_otu_paths(src)) == set(get_subdir(src))


@pytest.mark.parametrize("src", [SRC_PATH])
def test_utils_get_isolate_paths(src):
    assert set(get_isolate_paths(src)) == set(get_subdir(src))


def test_utils_hashing():
    initial_hashes = generate_unique_ids(n=40, length=8, mixed_case=False, excluded=[])

    additional_hashes = generate_unique_ids(
        n=20, length=8, mixed_case=False, excluded=initial_hashes
    )

    assert initial_hashes.isdisjoint(additional_hashes)
