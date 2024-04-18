import shutil
from io import StringIO
from pathlib import Path

import pytest
from rich.console import Console

from virtool_cli.legacy.utils import build_legacy_otu
from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ref.init import init_reference
from virtool_cli.ref.legacy import Repo
from virtool_cli.ref.repo import EventSourcedRepo
from virtool_cli.ref.utils import DataType


@pytest.fixture()
def console(mocker) -> Console:
    console = Console(file=StringIO())
    mocker.patch("virtool_cli.utils.console.console", console)
    return console


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

    return path


@pytest.fixture()
def scratch_repo(scratch_path: Path) -> Repo:
    """The path to a scratch legacy repository.

    TODO: Remove this along with flattened repo code.
    """
    return Repo(scratch_path)


@pytest.fixture()
def scratch_src_path(scratch_path: Path) -> Path:
    """The source path of a scratch reference repository.

    TODO: Remove this along with flattened repo code.
    """
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


@pytest.fixture()
def empty_repo(tmp_path: Path) -> EventSourcedRepo:
    return EventSourcedRepo.new(
        DataType.GENOME,
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


@pytest.fixture()
def legacy_repo_path(
    test_files_path: Path,
    tmp_path: Path,
) -> Path:
    """The path to a scratch legacy reference."""
    path = tmp_path / "reference"

    shutil.copytree(test_files_path / "src_v1", path / "src", dirs_exist_ok=True)

    return path


@pytest.fixture()
def legacy_otu(legacy_repo_path: Path) -> dict:
    return build_legacy_otu(
        legacy_repo_path / "src" / "a" / "abaca_bunchy_top_virus",
    )
