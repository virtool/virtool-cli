import shutil
from pathlib import Path

import pytest

from virtool_cli.ref.init import init_reference
from virtool_cli.ncbi.client import NCBIClient


@pytest.fixture()
def test_files_path():
    return Path(__file__).parent / "files"


@pytest.fixture()
def src_malformed_path(test_files_path: Path) -> Path:
    return test_files_path / "src_malformed"


@pytest.fixture()
def src_test_path(test_files_path: Path) -> Path:
    return test_files_path / "src_test"


@pytest.fixture()
def src_scratch_path(src_test_path: Path, tmp_path: Path) -> Path:
    path = tmp_path / "src_scratch"
    shutil.copytree(src_test_path, path)

    yield path

    shutil.rmtree(path)


@pytest.fixture()
def scratch_path(src_test_path: Path, cache_example_path: Path, tmp_path: Path) -> Path:
    path = tmp_path / "reference"
    init_reference(path)

    shutil.copytree(src_test_path, path / "src", dirs_exist_ok=True)
    shutil.copytree(cache_example_path, path / ".cache", dirs_exist_ok=True)

    yield path

    shutil.rmtree(path)


@pytest.fixture()
def cache_example_path(test_files_path: Path) -> Path:
    return test_files_path / "cache_test"


@pytest.fixture()
def cache_scratch_path(cache_example_path, tmp_path: Path) -> Path:
    path = tmp_path / "cache_scratch"
    shutil.copytree(cache_example_path, path)

    yield path

    shutil.rmtree(path)


@pytest.fixture()
def scratch_client(cache_scratch_path):
    return NCBIClient(cache_scratch_path, ignore_cache=False)
