import pytest
from pathlib import Path

from virtool_cli.repo.cls import Repo
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.model import NCBINuccore, NCBISource


ACCESSION_LISTS = [
    ["AB017504", "MH200607", "MK431779", "NC_003355"],
    ["NC_036587", "MT240513", "MT240490"],
]


@pytest.fixture()
def test_client(cache_scratch_path):
    return NCBIClient(cache_scratch_path / "cache_test")


@pytest.fixture
def test_records_path(test_files_path: Path):
    return test_files_path / "cache_test" / "nuccore"


@pytest.fixture
def test_taxonomy_path(test_files_path: Path):
    return test_files_path / "cache_test" / "taxonomy"


@pytest.fixture()
def scratch_repo(scratch_path):
    return Repo(scratch_path)


class TestClientFetchAccessions:
    @pytest.mark.parametrize("accessions", ACCESSION_LISTS)
    async def test_fetch_accessions_from_ncbi(self, accessions, test_client):
        clean_records = await test_client.fetch_accessions(
            accessions=accessions, use_cached=False
        )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.accession) is str

            assert type(record.source) is NCBISource

            assert type(record.source.taxid) is int

    @pytest.mark.parametrize("accessions", ACCESSION_LISTS)
    async def test_fetch_accessions_from_cache(self, accessions, test_client):
        clean_records = await test_client.fetch_accessions(
            accessions=accessions, use_cached=True
        )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.accession) is str

            assert type(record.source) is NCBISource

            assert type(record.source.taxid) is int

    @pytest.mark.asyncio
    async def test_fetch_accessions_fail(self, test_client):
        false_accessions = ["friday", "paella", "111"]

        records = await test_client.fetch_accessions(false_accessions)

        assert not records


@pytest.mark.parametrize("accessions", ACCESSION_LISTS)
class TestClientFetchRawAccessions:
    @pytest.mark.asyncio
    async def test_fetch_raw_via_accessions(self, accessions):
        records = await NCBIClient.fetch_raw_via_accessions(accessions)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_raw_via_accessions_partial(self, accessions):
        partial_accession_list = accessions
        partial_accession_list[0] = partial_accession_list[0][:3]

        records = await NCBIClient.fetch_raw_via_accessions(accessions)

        assert len(records) == len(accessions) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)


class TestClientTaxonomyUtilities:
    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxon_id", [908125, 1016856])
    async def test_fetch_taxonomy(self, taxon_id):
        taxonomy = await NCBIClient.fetch_taxonomy_by_taxid(taxon_id)

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
