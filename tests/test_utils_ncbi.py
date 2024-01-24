from pathlib import Path

import pytest

from virtool_cli.utils.ncbi import request_linked_accessions, fetch_taxonomy_rank


@pytest.fixture()
def taxon_ids(test_files_path: Path) -> list[int]:
    return [
        int(listing.stem.split("--")[0])
        for listing in (test_files_path / "catalog").glob("*--*.json")
    ]


@pytest.mark.parametrize("index", [0, 1, 2])
async def test_utils_request_linked_accessions(index: int, taxon_ids: list[int]):
    linked_accessions = await request_linked_accessions(taxon_ids[index])
    assert type(linked_accessions) is list

    for accession in linked_accessions:
        assert type(accession) == str


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
