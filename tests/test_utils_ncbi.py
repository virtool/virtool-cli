import pytest
from pathlib import Path

from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.utils.ncbi import request_linked_accessions, fetch_taxonomy_rank
from paths import TEST_FILES_PATH

CATALOG_PATH = TEST_FILES_PATH / "catalog"
SRC_PATH = TEST_FILES_PATH / "src_test"

TAXON_IDS = [listing.stem.split("--")[0] for listing in CATALOG_PATH.glob("*--*.json")]


@pytest.mark.asyncio
@pytest.mark.parametrize("index", [0, 1, 2])
async def test_utils_request_linked_accessions(index):
    linked_accessions = await request_linked_accessions(TAXON_IDS[index])
    assert type(linked_accessions) is list

    for accession in linked_accessions:
        assert type(accession) == str


@pytest.mark.asyncio
@pytest.mark.parametrize("index", [0, 1, 2])
async def test_utils_fetch_taxonomy_rank(index):
    taxon_id = TAXON_IDS[index]
    print(taxon_id)

    if taxon_id is "none":
        return

    rank = await fetch_taxonomy_rank(taxon_id)

    assert rank in ["no rank", "species", "genus", "subfamily", "family"]
