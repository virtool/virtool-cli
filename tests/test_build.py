import os
import json
import arrow
import pytest
import subprocess

from paths import TEST_FILES_PATH

TEST_SRC_PATH = TEST_FILES_PATH / "src_test"


@pytest.fixture()
def output_path(tmp_path):
    return tmp_path / "reference.json"


@pytest.fixture()
def command(output_path):
    return [
        "virtool",
        "ref",
        "build",
        "--output",
        str(output_path),
        "--src_path",
        str(TEST_SRC_PATH),
    ]


@pytest.mark.parametrize("version", [None, "v1.0.0", "v0.9.3"])
def test_version(version, command, output_path):
    """
    Test that the version field is correctly set in the reference.json file.
    """
    if version:
        command += ["--version", version]

    print(command)

    subprocess.run(command)
    with open(output_path, "r") as f:
        built_json = json.load(f)

    assert built_json["name"] == version


def test_created_at(command, output_path):
    """
    Test that the time of the creation in the reference.json file is correct
    """
    subprocess.run(command)
    with open(output_path, "r") as f:
        built_json = json.load(f)

    created_at = arrow.get(built_json["created_at"])

    assert (arrow.utcnow() - created_at).seconds == 0


# @pytest.mark.parametrize("indent", (True, False))
def test_indent(command, tmp_path):
    """
    Test that the indent in the reference.json file is properly set
    """
    output_path = tmp_path / "reference.py"
    output_indented_path = tmp_path / "reference_indent.py"

    subprocess.run([
        "virtool",
        "ref",
        "build",
        "--output",
        str(output_path),
        "--src_path",
        str(TEST_SRC_PATH),
    ])

    subprocess.run([
        "virtool",
        "ref",
        "build",
        "--output",
        str(output_indented_path),
        "--src_path",
        str(TEST_SRC_PATH),
        "--indent"
    ])

    expected_size = output_path.stat().st_size
    output_size = output_indented_path.stat().st_size

    assert expected_size != output_size
