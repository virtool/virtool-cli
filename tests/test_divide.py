import pytest
import subprocess
import filecmp

TEST_PATH = "tests/files/reference.json"
TEST_WITH_INDENT_PATH = "tests/files/reference_with_indent.json"
TEST_DIRECTORY_PATH = "tests/files/src"


@pytest.fixture()
def output(tmpdir):
    return tmpdir.join("src")


@pytest.fixture()
def command(output):
    return [
        "python", "virtool_cli/run.py",
        "divide", "-o", str(output),
        "-src"]


@pytest.mark.parametrize("src", [TEST_PATH,
                                 TEST_WITH_INDENT_PATH])
def test_divide(command, src, output, tmpdir):
    """
    Tests the divide operation to see if it produces a valid directory

    """
    command.append(src)

    subprocess.call(command)

    comparison = filecmp.dircmp(TEST_DIRECTORY_PATH, output)

    assert len(comparison.diff_files) == 0

