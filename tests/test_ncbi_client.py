import pytest
from pathlib import Path

from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.utils import NuccorePacket
from virtool_cli.ncbi.model import NCBIAccession, NCBISource, NCBITaxonomy, NCBIRank


@pytest.fixture()
def test_client(test_files_path):
    return NCBIClient(test_files_path / "cache_test")


@pytest.fixture
def test_records_path(test_files_path: Path):
    return test_files_path / "cache_test" / "nuccore"


@pytest.fixture
def test_taxonomy_path(test_files_path: Path):
    return test_files_path / "cache_test" / "taxonomy"


@pytest.fixture()
def accession_list():
    return ["KT390494", "KT390496", "KT390501"]


class TestClientMain:
    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxon_id", [1016856])
    async def test_procure_from_taxonomy(self, taxon_id, test_client):
        clean_records = await test_client.procure_from_taxonomy(
            taxon_id, use_cached=False
        )

        assert type(clean_records) is list

        for packet in clean_records:
            assert type(packet) is NuccorePacket

            assert type(packet.sequence) is NCBIAccession

            assert type(packet.source) is NCBISource

            assert type(packet.sequence.accession) is str

            assert type(packet.source.taxid) is int


class TestClientUtilities:
    @pytest.mark.asyncio
    async def test_fetch_accessions(self, accession_list):
        records = await NCBIClient.fetch_accessions(accession_list)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_accessions_partial(self, accession_list):
        partial_accession_list = accession_list
        partial_accession_list[0] = partial_accession_list[0][:3]

        records = await NCBIClient.fetch_accessions(accession_list)

        assert len(records) == len(accession_list) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_accessions_fail(self):
        accession_list = ["friday", "paella", "111"]

        assert not await NCBIClient.fetch_accessions(accession_list)

    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxon_id", [908125, 1016856])
    async def test_fetch_taxonomy(self, taxon_id):
        taxonomy = await NCBIClient.fetch_taxonomy(taxon_id)

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
    async def test_fetch_taxonomy_by_name(self, name, taxid):
        assert await NCBIClient.fetch_taxonomy_id_by_name(name) == taxid

    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "taxid,rank",
        [
            (438782, "species"),
            (1016856, "isolate"),
        ],
    )
    async def test_fetch_taxonomy_rank(self, taxid, rank):
        fetched_rank = await NCBIClient.fetch_taxon_rank(taxid)

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
    async def test_check_spelling(self, misspelled, expected):
        taxon_name = await NCBIClient.check_spelling(name=misspelled)

        assert taxon_name.lower() == expected.lower()
