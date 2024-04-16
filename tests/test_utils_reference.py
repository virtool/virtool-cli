from pathlib import Path

import pytest

from virtool_cli.utils.id_generator import generate_unique_ids
from virtool_cli.utils.reference import (
    get_isolate_paths,
    get_otu_paths,
    is_v1,
    search_otu_by_id,
)


def get_subdir(path: Path) -> list:
    """Generates a list of all subdirectories under a path

    :param path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    return [subdir for subdir in path.iterdir() if subdir.is_dir()]


def test_utils_is_v1(test_files_path: Path):
    assert is_v1(test_files_path / "src_v1")


def test_utils_is_v2(scratch_src_path: Path):
    assert not is_v1(scratch_src_path)


@pytest.mark.parametrize("otu_id", ["f8a56910", "3962f6ec"])
def test_utils_search_by_id_success(otu_id: str, scratch_src_path: Path):
    assert search_otu_by_id(otu_id, scratch_src_path) is not None


@pytest.mark.parametrize("otu_id", ["bbbbbbbb", "777777tt"])
def test_utils_search_by_id_fail(otu_id: str, scratch_src_path: Path):
    assert search_otu_by_id(otu_id, scratch_src_path) is None


def test_utils_get_otu_paths(scratch_src_path: Path):
    assert set(get_otu_paths(scratch_src_path)) == set(get_subdir(scratch_src_path))


def test_utils_get_isolate_paths(scratch_src_path: Path):
    assert set(get_isolate_paths(scratch_src_path)) == set(get_subdir(scratch_src_path))


def test_utils_hashing():
    initial_hashes = generate_unique_ids(n=40, length=8, mixed_case=False, excluded=[])

    additional_hashes = generate_unique_ids(
        n=20,
        length=8,
        mixed_case=False,
        excluded=initial_hashes,
    )

    assert initial_hashes.isdisjoint(additional_hashes)
