import shutil
import subprocess
from pathlib import Path

import pytest


def test_empty_success(tmp_path: Path):
    subprocess.run(["virtool", "ref", "init", "-repo", str(tmp_path)])

    p = subprocess.run(
        [
            "virtool",
            "ref",
            "update",
            "reference",
            "-src",
            str(tmp_path / "src"),
        ],
        capture_output=True,
    )

    assert p.returncode == 0


@pytest.mark.parametrize(
    "otu_dirname",
    [
        "abaca_bunchy_top_virus--c93ec9a9",
        "babaco_mosaic_virus--xcl20vqt",
        "cabbage_leaf_curl_jamaica_virus--d226290f",
        "faba_bean_necrotic_stunt_alphasatellite_1--6444acf3",
    ],
)
def test_update_otu(otu_dirname: str, src_test_path: Path, tmp_path: Path):
    src_path = tmp_path / "src"
    shutil.copytree(src_test_path, src_path)

    otu_path = src_path / otu_dirname

    subprocess.run(
        [
            "virtool",
            "ref",
            "update",
            "otu",
            "-otu",
            str(otu_path),
        ]
    )

    assert set(otu_path.iterdir()) - set((src_test_path / otu_dirname).iterdir())
