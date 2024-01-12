import pytest

from virtool_cli.utils.storage import get_otu_accessions
from paths import TEST_FILES_PATH

SRC_PATH = TEST_FILES_PATH / "src_test"

test_dirs = [
    "abaca_bunchy_top_virus--c93ec9a9",
    "babaco_mosaic_virus--xcl20vqt",
    "cabbage_leaf_curl_jamaica_virus--d226290f",
    "faba_bean_necrotic_stunt_alphasatellite_1--6444acf3",
]

@pytest.mark.asyncio
@pytest.mark.parametrize("otu_dirname", test_dirs)
async def test_utils_get_otu_accessions(otu_dirname):
    accession_list = await get_otu_accessions(SRC_PATH / otu_dirname)
    assert accession_list != []
