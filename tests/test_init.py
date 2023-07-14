import pytest
import subprocess
import filecmp

from paths import TEST_FILES_PATH

TEST_PATH = TEST_FILES_PATH / "repo_init"
TEST_COMP_PATH = TEST_FILES_PATH / "repo_ex"

@pytest.fixture()
def command(output):
    return ["python", "virtool_cli/run.py", "ref", "init", "-repo", TEST_PATH]

@pytest.mark.parametrize("src", [TEST_PATH])
def test_init(command, src, output, tmpdir):
    """
    Tests the divide init to see if it produces a valid directory

    """
    command.append(src)

    subprocess.call(command)

    comparison = filecmp.dircmp(TEST_COMP_PATH, output)

    assert len(comparison.diff_files) == 0

    subprocess.call(["rm", "-r", TEST_PATH])

