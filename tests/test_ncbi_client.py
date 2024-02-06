import json
from pathlib import Path

import pytest

from virtool_cli.ncbi.client import NCBIClient


@pytest.fixture()
def client():
    return NCBIClient()


class TestClient:
    @pytest.fixture()
    def accession_list(self):
        return ["KT390494", "KT390496", "KT390501"]

    @pytest.mark.asyncio
    async def test_fetch_accessions(self, client, accession_list):
        records = await client.fetch_accessions(accession_list)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxon_id", [908125, 1016856])
    async def test_fetch_taxonomy(self, client, taxon_id):
        taxonomy = await client.fetch_taxonomy(taxon_id, long=True)

        assert taxonomy.get("TaxId", None) is not None
        assert taxonomy.get("ScientificName", None) is not None
        assert taxonomy.get("Rank", None) is not None

    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "name,taxid",
        [
            ("Abaca bunchy top virus", 438782),
            ("Rhynchosia golden mosaic virus", 117198),
        ],
    )
    async def test_fetch_taxonomy_by_name(self, client, name, taxid):
        taxon_id_by_name = await client.fetch_taxonomy_id_by_name(name)

        assert taxon_id_by_name == taxid

    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "taxid,rank",
        [
            (438782, "species"),
            (1016856, "isolate"),
        ],
    )
    async def test_fetch_taxonomy_rank(self, client, taxid, rank):
        fetched_rank = await client.fetch_taxon_rank(taxid)

        assert fetched_rank == rank

    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "misspelled,expected",
        [
            (
                "Hynchosia yellow mosaic India virus",
                "Rhynchosia yellow mosaic India virus",
            ),
            ("Angelica bush stunt virus", "Angelica bushy stunt virus"),
        ],
    )
    async def test_check_spelling(self, client, misspelled, expected):
        taxon_name = await client.check_spelling(name=misspelled)

        assert taxon_name.lower() == expected.lower()
