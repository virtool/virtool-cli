import pytest
import subprocess
from paths import TEST_FILES_PATH

TEST_SRC_PATH = TEST_FILES_PATH / "src_test"
TEST_BAD_SRC_PATH = TEST_FILES_PATH / "src_malformed"


@pytest.fixture()
def command(otu_path):
    return ["virtool", "ref", "check", "otu", "--otu_path", str(otu_path)]


@pytest.mark.parametrize("otu_dir", ["abaca_bunchy_top_virus--c93ec9a9"])
def test_otu_check_success(otu_dir):
    otu_path = TEST_SRC_PATH / otu_dir
    output = subprocess.check_output(
        ["virtool", "ref", "check", "otu", "--otu_path", str(otu_path)]
    ).decode("utf-8")

    assert output == ""


@pytest.mark.parametrize("otu_dir", ["abaca_bunchy_top_virus--c93ec9a9"])
def test_otu_check_fail(otu_dir):
    otu_path = TEST_BAD_SRC_PATH / otu_dir
    output = subprocess.check_output(
        ["virtool", "ref", "check", "otu", "--otu_path", str(otu_path)]
    ).decode("utf-8")

    assert output != ""
