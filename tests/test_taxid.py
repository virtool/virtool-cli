import pytest
import json
import subprocess
import shutil

from paths import TEST_FILES_PATH
TEST_PATH = TEST_FILES_PATH / "src_notaxid"

# @pytest.mark.skip(reason="Moved to accession catalog")
@pytest.mark.parametrize("base_path", [TEST_PATH])
def test_fetch_taxids(base_path, tmp_path):
    """
    """
    fetch_path = tmp_path / 'src'

    shutil.copytree(base_path, fetch_path)
    assert len([fetch_path.glob('[a-z]')]) > 0

    subprocess.call([
        "virtool", "ref", "doctor", "taxid",
        "-src", str(fetch_path),
        "--force_update"
    ])

    for otu_path in fetch_path.iterdir():
        if not otu_path.is_dir():
            continue

        with open(otu_path / "otu.json", "r") as f:
            otu = json.load(f)
        
        assert otu['taxid'] is not None