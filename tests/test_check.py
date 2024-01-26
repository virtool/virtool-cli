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
def malformed_src(src_scratch_path):
    # Malformation 1: Remove all isolates
    for isolate_path in get_isolate_paths(src_scratch_path / OTU_DIR_NAMES[0]):
        shutil.rmtree(isolate_path)

    # Malformation 2: Remove isolate metadata file
    for isolate_path in get_isolate_paths(src_scratch_path / OTU_DIR_NAMES[1]):
        (isolate_path / "isolate.json").unlink()

    # Malformation 3: Remove critical value from OTU metadata
    otu_path = src_scratch_path / OTU_DIR_NAMES[2]
    otu = json.loads((otu_path / 'otu.json').read_text())
    otu["_id"] = ""
    with open(otu_path / 'otu.json', "w") as f:
        json.dump(otu, f, indent=4)

    # Malformation 4: Remove critical value from isolate metadata
    for isolate_path in get_isolate_paths(src_scratch_path / OTU_DIR_NAMES[3]):
        isolate = json.loads((isolate_path / 'isolate.json').read_text())
        isolate["id"] = ""
        with open(isolate_path / 'isolate.json', "w") as f:
            json.dump(isolate, f, indent=4)

    return src_scratch_path


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
def test_otu_check_fail(otu_dir: str, malformed_src: Path):
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
