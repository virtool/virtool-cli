import os
import json
import arrow
import pytest
import subprocess

TEST_PATH = "tests/files/reference.json"
TEST_WITH_INDENT_PATH = "tests/files/reference_with_indent.json"


@pytest.fixture()
def output(tmpdir):
    return tmpdir.join("reference.json")


@pytest.fixture()
def command(output):
    return [
        "virtool", "build", "-o", str(output),
        "-src", "tests/files/src"]


@pytest.mark.parametrize("version", [None, "v1.0.0", "v0.9.3"])
def test_version(version, command, output):
    """
    Test that the version field is correctly set in the reference.json file.

    """

    if version:
        command += ["-V", version]

    subprocess.call(command)

    built_json = json.load(output)

    assert built_json["name"] == version


def test_created_at(command, output):
    """
    Test that the time of the creation in the reference.json file is correct

    """

    subprocess.call(command)

    built_json = json.load(output)

    created_at = arrow.get(built_json["created_at"])

    assert (arrow.utcnow() - created_at).seconds == 0


@pytest.mark.parametrize("indent", (True, False))
def test_indent(command, output, indent):
    """
    Test that the indent in the reference.json file is properly set

    """
    generate_expected = [
        "virtool", "build", "-src",
        "tests/files/src"]
    if indent:
        command.append("-i")
        generate_expected.append("-i")
        expected_path = TEST_WITH_INDENT_PATH
    else:
        expected_path = TEST_PATH

    generate_expected.extend(["-o", expected_path])

    subprocess.call(generate_expected)
    subprocess.call(command)

    expected_size = os.path.getsize(expected_path)
    output_size = os.path.getsize(output)
    print(expected_size, output_size)

    assert abs(expected_size - output_size) < 10
