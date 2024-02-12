import pytest
from urllib.error import HTTPError

from virtool_cli.ncbi.client import NCBIClient


@pytest.fixture()
def test_client():
    return NCBIClient()


@pytest.fixture()
def accession_list():
    return ["KT390494", "KT390496", "KT390501"]


class TestClient:
    @pytest.mark.asyncio
    async def test_fetch_accessions(self, test_client, accession_list):
        records = await test_client.fetch_accessions(accession_list)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_accessions_partial(self, test_client, accession_list):
        partial_accession_list = accession_list
        partial_accession_list[0] = partial_accession_list[0][:3]

        records = await test_client.fetch_accessions(accession_list)

        assert len(records) == len(accession_list) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_accessions_fail(self, test_client):
        accession_list = ["friday", "paella", "111"]

        with pytest.raises(HTTPError):
            await test_client.fetch_accessions(accession_list)

    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxon_id", [908125, 1016856])
    async def test_fetch_taxonomy(self, test_client, taxon_id):
        taxonomy = await test_client.fetch_taxonomy(taxon_id, long=True)

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
    async def test_fetch_taxonomy_by_name(self, test_client, name, taxid):
        assert await test_client.fetch_taxonomy_id_by_name(name) == taxid

    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "taxid,rank",
        [
            (438782, "species"),
            (1016856, "isolate"),
        ],
    )
    async def test_fetch_taxonomy_rank(self, test_client, taxid, rank):
        fetched_rank = await test_client.fetch_taxon_rank(taxid)

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
    async def test_check_spelling(self, test_client, misspelled, expected):
        taxon_name = await test_client.check_spelling(name=misspelled)

        assert taxon_name.lower() == expected.lower()
