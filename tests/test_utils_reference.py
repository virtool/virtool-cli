import pytest

from virtool_cli.utils.reference import *
from paths import TEST_FILES_PATH
SRC_PATH = TEST_FILES_PATH / "src_test"
SRC_V1_PATH = TEST_FILES_PATH / "src_v1"


def test_utils_is_v1():
    assert is_v1(SRC_V1_PATH)

def test_utils_is_v2():
    assert not is_v1(SRC_PATH)

@pytest.mark.parametrize("otu_id", ['f8a56910', '3962f6ec'])
def test_utils_search_by_id(otu_id):
    assert search_otu_by_id(SRC_PATH, otu_id) is not None