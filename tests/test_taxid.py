import pytest
import os
import json
import subprocess

TEST_PATH = "tests/files/src_a"


@pytest.fixture(scope="session", autouse=True)
def command():
    command = [
        "python", "virtool_cli/run.py",
        "taxid", "-src",
        TEST_PATH, "-f"]
    subprocess.call(command)


@pytest.mark.parametrize("path", ["h/hop_stunt_viroid", "r/reovirus_tf1_(not_a_plant_virus)",
                                  "t/tobacco_mosaic_virus", "t/totivirus_tf1_(not_a_plant_virus)"])
def test_taxid(path, command):
    path = os.path.join(TEST_PATH, path)

    with open(os.path.join(path, "otu.json"), 'r') as f:
        otu = json.load(f)

        if "not" in otu["name"]:
            assert otu["taxid"] is None
        else:
            assert otu["taxid"] is not None
