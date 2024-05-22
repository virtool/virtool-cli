import json
import subprocess
from pathlib import Path

import arrow
import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from virtool_cli.ref.build import build_json


@pytest.fixture()
def output_path(tmp_path):
    return tmp_path / "reference.json"


def test_ok(output_path, scratch_path: Path, snapshot: SnapshotAssertion):
    """Test that the command exits with a 0 exit code."""
    build_json(False, output_path, scratch_path, "2.1.0")

    with open(output_path) as f:
        built_json = json.load(f)

    assert built_json == snapshot(exclude=props("created_at"))
    assert built_json["name"] == "2.1.0"
    assert (arrow.utcnow() - arrow.get(built_json["created_at"])).seconds == 0


def test_indent(scratch_path: Path, tmp_path: Path):
    """Test that the indent in the reference.json file is properly set"""
    output_path = tmp_path / "reference.json"
    output_indented_path = tmp_path / "reference_indent.json"

    build_json(False, output_path, scratch_path, "2.1.0")
    build_json(True, output_indented_path, scratch_path, "2.1.0")

    assert {**json.load(open(output_path, "rb")), "created_at": ""} == {
        **json.load(open(output_indented_path, "rb")),
        "created_at": "",
    }
