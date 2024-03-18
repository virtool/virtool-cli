import pytest
from syrupy import SnapshotAssertion
from structlog import get_logger

from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.ncbi.model import NCBINuccore, NCBISource, NCBITaxonomy

test_logger = get_logger()

SAFETY_PAUSE = 3


@pytest.fixture()
def scratch_client(cache_scratch_path):
    return NCBIClient(cache_scratch_path)


@pytest.fixture()
def empty_client(tmp_path):
    return NCBIClient(NCBICache(tmp_path / "clean_repo").path)


class TestClientFetchGenbank:
    @pytest.mark.ncbi
    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    async def test_fetch_genbank_records_from_ncbi(
        self, accessions, empty_client, snapshot: SnapshotAssertion
    ):
        clean_records = await empty_client.fetch_genbank_records(
            accessions=accessions, cache_results=True, use_cached=False
        )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.source) is NCBISource

        assert {
            path.name for path in empty_client.cache._nuccore_path.glob("*.json")
        } == snapshot

    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    async def test_fetch_genbank_records_from_cache(
        self, accessions, cache_scratch_path, snapshot: SnapshotAssertion
    ):
        client = NCBIClient(cache_scratch_path)

        clean_records = await client.fetch_genbank_records(
            accessions=accessions, cache_results=False, use_cached=True
        )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.source) is NCBISource

            assert record == snapshot(name=f"{record.accession}_validated.json")

    @pytest.mark.ncbi
    @pytest.mark.parametrize(
        "accessions", [["AB017503", "AB017504", "MH200607", "MK431779", "NC_003355"]]
    )
    async def test_fetch_partially_cached_genbank_records(
        self, accessions: list[str], cache_scratch_path
    ):
        client = NCBIClient(cache_scratch_path)

        try:
            assert next(client.cache._nuccore_path.glob("*.json"))
        except StopIteration:
            pytest.fail("Could not retrieve any records from cache")

        clean_records = await client.fetch_genbank_records(
            accessions=accessions, cache_results=False, use_cached=True
        )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.source) is NCBISource

    @pytest.mark.ncbi
    @pytest.mark.asyncio
    async def test_fetch_accessions_fail(self, scratch_client):
        false_accessions = ["friday", "paella", "111"]

        records = await scratch_client.fetch_genbank_records(false_accessions)

        assert not records


class TestClientFetchRawGenbank:
    @pytest.mark.ncbi
    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    async def test_fetch_raw_via_accessions(self, accessions: list[str]):
        records = await NCBIClient.fetch_unvalidated_genbank_records(accessions)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.ncbi
    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "accessions",
        [
            ["paella", "MH200607", "MK431779", "NC_003355"],
            ["friday", "MT240513", "MT240490"],
        ],
    )
    async def test_fetch_raw_via_accessions_partial(self, accessions: list[str]):
        records = await NCBIClient.fetch_unvalidated_genbank_records(accessions)

        assert len(records) == len(accessions) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.ncbi
    @pytest.mark.asyncio
    async def test_fetch_raw_via_accessions_fail(self):
        false_accessions = ["friday", "paella", "111"]
        records = await NCBIClient.fetch_unvalidated_genbank_records(false_accessions)

        assert not records


@pytest.mark.ncbi
@pytest.mark.asyncio
@pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
async def test_fetch_records_by_taxid(taxid, empty_client, snapshot):
    with pytest.raises(StopIteration):
        next(empty_client.cache._nuccore_path.glob("*.json"))

    records = await empty_client.link_from_taxid_and_fetch(taxid, cache_results=True)

    assert records

    for record in records:
        assert type(record) is NCBINuccore

    assert {
        path.name for path in empty_client.cache._nuccore_path.glob("*.json")
    } == snapshot


class TestClientFetchTaxonomy:
    @pytest.mark.ncbi
    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
    async def test_fetch_taxonomy_from_ncbi(
        self, taxid, empty_client, snapshot: SnapshotAssertion
    ):
        taxonomy = await empty_client.fetch_taxonomy_record(
            taxid, use_cached=False, cache_results=True
        )

        assert type(taxonomy) is NCBITaxonomy

        assert {
            path.name for path in empty_client.cache._taxonomy_path.glob("*.json")
        } == snapshot

    @pytest.mark.ncbi
    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxid", [1000000000000, 99999999])
    async def test_fetch_taxonomy_from_ncbi_fail(self, taxid: int, empty_client):
        taxonomy = await empty_client.fetch_taxonomy_record(taxid, use_cached=False)

        assert taxonomy is None

    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxid", [438782, 1198450])
    async def test_fetch_taxonomy_from_cache(
        self, taxid: int, scratch_client, snapshot: SnapshotAssertion
    ):
        assert scratch_client.cache.load_taxonomy(taxid)

        taxonomy = await scratch_client.fetch_taxonomy_record(taxid, use_cached=True)

        assert type(taxonomy) is NCBITaxonomy

        assert taxonomy == snapshot(name=f"{taxonomy}_validated.json")


@pytest.mark.ncbi
@pytest.mark.asyncio
@pytest.mark.parametrize(
    "name,taxid",
    [
        ("Abaca bunchy top virus", 438782),
        ("Rhynchosia golden mosaic virus", 117198),
    ],
)
async def test_fetch_taxonomy_by_name(name: str, taxid: int):
    assert await NCBIClient.fetch_taxonomy_id_by_name(name) == taxid


@pytest.mark.ncbi
@pytest.mark.asyncio
@pytest.mark.parametrize(
    "misspelled,expected",
    [
        (
            "Hynchosia yellow mosaic India virus",
            "Rhynchosia yellow mosaic India virus",
        ),
        (
            "Angelica bush stunt virus",
            "Angelica bushy stunt virus",
        ),
    ],
)
async def test_check_spelling(misspelled: str, expected: str):
    taxon_name = await NCBIClient.check_spelling(name=misspelled)

    assert taxon_name.lower() == expected.lower()
