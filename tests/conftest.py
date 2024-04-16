import shutil
from pathlib import Path

import pytest

from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ref.init import init_reference


@pytest.fixture()
def test_files_path():
    return Path(__file__).parent / "files"


@pytest.fixture()
def scratch_path(
    test_files_path: Path,
    tmp_path: Path,
) -> Path:
    """The path to a scratch reference repository."""
    path = tmp_path / "reference"

    init_reference(path)

    shutil.copytree(test_files_path / "src_test", path / "src", dirs_exist_ok=True)
    shutil.copytree(test_files_path / "cache_test", path / ".cache", dirs_exist_ok=True)

    yield path

    shutil.rmtree(path)


@pytest.fixture()
def scratch_src_path(scratch_path: Path) -> Path:
    """The source path of a scratch reference repository."""
    return scratch_path / "src"


@pytest.fixture()
def scratch_ncbi_cache_path(
    scratch_path: Path,
) -> Path:
    """The path to a scratch NCBI client cache."""
    return scratch_path / ".cache"


@pytest.fixture()
def scratch_ncbi_cache(scratch_ncbi_cache_path: Path):
    """A scratch NCBI cache with preloaded data."""
    return NCBICache(scratch_ncbi_cache_path)


@pytest.fixture()
def scratch_ncbi_client(scratch_ncbi_cache_path: Path):
    """A scratch NCBI client with a preloaded cache."""
    return NCBIClient(scratch_ncbi_cache_path, ignore_cache=False)
