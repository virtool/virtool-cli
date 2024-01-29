import pytest

from virtool_cli.utils.ncbi import request_linked_accessions, fetch_taxonomy_rank
from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.utils.storage import read_otu


@pytest.fixture()
async def taxon_ids(src_test_path) -> list[int]:
    taxon_ids = []

    for otu_path in get_otu_paths(src_test_path):
        otu = await read_otu(otu_path)
        taxon_ids.append(otu["taxid"])

    return taxon_ids


@pytest.mark.parametrize("index", [0, 1, 2])
async def test_utils_request_linked_accessions(index: int, taxon_ids: list[int]):
    accessions = await request_linked_accessions(taxon_ids[index])

    assert isinstance(accessions, list)
    assert all(isinstance(accession, str) for accession in accessions)


@pytest.mark.skip
@pytest.mark.parametrize("index", [0, 1, 2])
async def test_utils_fetch_taxonomy_rank(index: int, taxon_ids: list[int]):
    taxon_id = taxon_ids[index]

    if taxon_id != "none":
        assert await fetch_taxonomy_rank(taxon_id) in [
            "no rank",
            "species",
            "genus",
            "subfamily",
            "family",
        ]
