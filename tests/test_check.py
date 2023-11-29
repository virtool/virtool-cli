import pytest
import subprocess
from paths import TEST_FILES_PATH

TEST_SRC_PATH = TEST_FILES_PATH / "src_test"
TEST_BAD_SRC_PATH = TEST_FILES_PATH / "src_malformed"

TEST_DIRS = [
    "abaca_bunchy_top_virus--c93ec9a9",
    "babaco_mosaic_virus--xcl20vqt",
    "cabbage_leaf_curl_jamaica_virus--d226290f",
    "faba_bean_necrotic_stunt_alphasatellite_1--6444acf3",
]


@pytest.fixture()
def command(otu_path):
    return ["virtool", "ref", "check", "otu", "--otu_path", str(otu_path)]


@pytest.mark.parametrize("otu_dir", TEST_DIRS)
def test_otu_check_success(otu_dir):
    otu_path = TEST_SRC_PATH / otu_dir
    output = subprocess.check_output(
        ["virtool", "ref", "check", "otu", "--otu_path", str(otu_path)]
    ).decode("utf-8")

    assert output.strip() == "True"


@pytest.mark.parametrize("otu_dir", TEST_DIRS)
def test_otu_check_fail(otu_dir):
    otu_path = TEST_BAD_SRC_PATH / otu_dir
    output = subprocess.check_output(
        ["virtool", "ref", "check", "otu", "--otu_path", str(otu_path)]
    ).decode("utf-8")

    assert output.strip() == "False"
