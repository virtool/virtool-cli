import json
from pathlib import Path

import pytest

from virtool_cli.ref.init import init_reference


@pytest.fixture()
def assert_repo_valid():
    """Provides a function that asserts that a reference repository is valid."""

    def func(path):
        assert (path / "src").exists()
        assert (path / ".cache").exists()
        assert (path / ".github").exists()

        with open(path / "src" / "meta.json") as f:
            assert json.load(f) == {"data_type": "genome"}

    return func


def test_init_empty(assert_repo_valid, log, tmp_path: Path):
    """Test that a new reference repository is created when the target path does not
    exist.
    """
    path = tmp_path / "new_repo"

    init_reference(path)

    assert_repo_valid(path)

    assert log.has("Complete")


def test_init_exists(assert_repo_valid, log, tmp_path: Path):
    """Test that a new reference repository can be created in an existing empty
    directory.
    """
    path = tmp_path / "new_repo"
    path.mkdir()

    init_reference(path)

    assert_repo_valid(path)

    assert log.has("Complete")


def test_init_not_empty(log, tmp_path: Path):
    """Test that reference initialization fails when the target path is not empty."""
    path = tmp_path / "new_repo"
    path.mkdir()
    (path / "file.txt").touch()

    init_reference(path)

    assert log.has("The directory is not empty")
    assert not log.has("Complete")


def test_init_is_file(log, tmp_path: Path):
    """Test that reference initialization fails when the target path is a file."""
    path = tmp_path / "new_repo"
    path.touch()

    init_reference(path)

    assert log.has("The target path is a file")
    assert not log.has("Complete")
