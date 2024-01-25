import shutil
import json
from pathlib import Path

import pytest
import subprocess

from virtool_cli.utils.reference import get_isolate_paths


OTU_DIR_NAMES = [
    "abaca_bunchy_top_virus--c93ec9a9",
    "babaco_mosaic_virus--xcl20vqt",
    "cabbage_leaf_curl_jamaica_virus--d226290f",
    "faba_bean_necrotic_stunt_alphasatellite_1--6444acf3",
]


@pytest.fixture()
def malformed_src(tmp_path, src_scratch_path):
    test_src_path = tmp_path / "src"
    shutil.copytree(src_scratch_path, test_src_path)

    delete_isolates(test_src_path / OTU_DIR_NAMES[0])

    remove_otu_id(test_src_path / OTU_DIR_NAMES[1])

    delete_isolate_metadata(test_src_path / OTU_DIR_NAMES[2])

    remove_isolate_id(test_src_path / OTU_DIR_NAMES[3])

    return test_src_path


@pytest.mark.parametrize("otu_dir", OTU_DIR_NAMES)
def test_otu_check_success(otu_dir: str, src_scratch_path: Path):
    output = subprocess.check_output(
        [
            "virtool",
            "ref",
            "check",
            "otu",
            "--otu_path",
            str(src_scratch_path / otu_dir),
        ]
    ).decode("utf-8")

    assert output.strip() == "True"


@pytest.mark.parametrize("otu_dir", OTU_DIR_NAMES)
def test_otu_check_fail(otu_dir: str, malformed_src: Path, tmp_path: Path):
    otu_path = malformed_src / otu_dir

    output = subprocess.check_output(
        [
            "virtool",
            "ref",
            "check",
            "otu",
            "--otu_path",
            str(otu_path),
        ]
    ).decode("utf-8")

    assert output.strip() == "False"


def delete_isolates(otu_path):
    for isolate_path in get_isolate_paths(otu_path):
        shutil.rmtree(isolate_path)


def delete_isolate_metadata(otu_path):
    for isolate_path in get_isolate_paths(otu_path):
        (isolate_path / "isolate.json").unlink()


def remove_otu_id(otu_path):
    otu = json.loads((otu_path / 'otu.json').read_text())
    otu["_id"] = ""
    with open(otu_path / 'otu.json', "w") as f:
        json.dump(otu, f, indent=4)


def remove_isolate_id(otu_path):
    for isolate_path in get_isolate_paths(otu_path):
        isolate = json.loads((isolate_path / 'isolate.json').read_text())
        isolate["id"] = ""
        with open(isolate_path / 'isolate.json', "w") as f:
            json.dump(isolate, f, indent=4)