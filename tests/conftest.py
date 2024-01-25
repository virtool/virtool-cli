import shutil
from pathlib import Path

import pytest

from virtool_cli.ref.init import init_reference


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def test_files_path():
    return Path(__file__).parent / "files"


@pytest.fixture()
def src_malformed_path(test_files_path: Path):
    return test_files_path / "src_malformed"


@pytest.fixture()
def src_test_path(test_files_path: Path):
    return test_files_path / "src_test"


@pytest.fixture()
def src_scratch_path(src_test_path: Path, tmp_path: Path):
    path = tmp_path / "src_scratch"
    shutil.copytree(src_test_path, path)

    return path


@pytest.fixture()
def scratch_path(src_test_path: Path, tmp_path: Path):
    path = tmp_path / "reference"
    init_reference(path)

    shutil.copytree(src_test_path, path / "src", dirs_exist_ok=True)

    return path
