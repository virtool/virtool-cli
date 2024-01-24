import shutil
from pathlib import Path

import pytest
import subprocess


OTU_DIR_NAMES = [
    "abaca_bunchy_top_virus--c93ec9a9",
    "babaco_mosaic_virus--xcl20vqt",
    "cabbage_leaf_curl_jamaica_virus--d226290f",
    "faba_bean_necrotic_stunt_alphasatellite_1--6444acf3",
]


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
def test_otu_check_fail(otu_dir: str, src_malformed_path: Path, tmp_path: Path):
    malformed_scratch_path = tmp_path / "malformed_scratch"
    shutil.copytree(src_malformed_path, malformed_scratch_path)

    output = subprocess.check_output(
        [
            "virtool",
            "ref",
            "check",
            "otu",
            "--otu_path",
            str(malformed_scratch_path / otu_dir),
        ]
    ).decode("utf-8")

    assert output.strip() == "False"
