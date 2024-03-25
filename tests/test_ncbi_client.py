from shutil import rmtree
from structlog import get_logger

import pytest
from syrupy import SnapshotAssertion

from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.model import NCBINuccore, NCBISource, NCBITaxonomy

test_logger = get_logger()

SERVER_ERRORS = ["HTTPError", "IncompleteRead"]


@pytest.fixture()
def scratch_client(cache_scratch_path):
    return NCBIClient(cache_scratch_path, ignore_cache=False)


@pytest.fixture()
def empty_client(tmp_path):
    dummy_cache_path = tmp_path / "dummy_cache"
    dummy_cache_path.mkdir()

    yield NCBIClient(dummy_cache_path, ignore_cache=True)

    rmtree(dummy_cache_path)


class TestClientFetchGenbank:
    @pytest.mark.ncbi
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
        clean_records = await empty_client.fetch_genbank_records(accessions=accessions)

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
        client = NCBIClient(cache_scratch_path, ignore_cache=False)

        clean_records = await client.fetch_genbank_records(accessions=accessions)

        assert clean_records

        for record in clean_records:
            assert record == snapshot(name=f"{record.accession}_validated.json")

    @pytest.mark.ncbi
    @pytest.mark.parametrize(
        "accessions", [["AB017503", "AB017504", "MH200607", "MK431779", "NC_003355"]]
    )
    async def test_fetch_partially_cached_genbank_records(
        self, accessions: list[str], cache_scratch_path, snapshot: SnapshotAssertion
    ):
        client = NCBIClient(cache_scratch_path, ignore_cache=False)

        assert list(client.cache._nuccore_path.glob("*.json")) != []

        clean_records = await client.fetch_genbank_records(accessions=accessions)

        assert clean_records

        for record in clean_records:
            assert record == snapshot(name=f"{record.accession}_validated.json")

    @pytest.mark.ncbi
    async def test_fetch_accessions_fail(self, scratch_client):
        records = await scratch_client.fetch_genbank_records(
            ["friday", "paella", "111"]
        )

        assert not records


class TestClientFetchRawGenbank:
    @pytest.mark.ncbi
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
            assert record.get("GBSeq_locus")
            assert record.get("GBSeq_sequence")

    @pytest.mark.ncbi
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
    async def test_fetch_raw_via_accessions_fail(self):
        records = await NCBIClient.fetch_unvalidated_genbank_records(
            ["friday", "paella", "111"]
        )

        assert not records


@pytest.mark.ncbi
@pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
async def test_fetch_records_by_taxid(taxid, empty_client, snapshot):
    assert not list(empty_client.cache._nuccore_path.glob("*.json"))

    records = await empty_client.link_from_taxid_and_fetch(taxid)

    assert records

    for record in records:
        assert type(record) is NCBINuccore

    assert {
        path.name for path in empty_client.cache._nuccore_path.glob("*.json")
    } == snapshot


class TestClientFetchTaxonomy:
    @pytest.mark.ncbi
    @pytest.mark.flaky(reruns=3, reruns_delay=3, only_rerun=SERVER_ERRORS)
    @pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
    async def test_fetch_taxonomy_from_ncbi(
        self, taxid, empty_client, snapshot: SnapshotAssertion
    ):
        taxonomy = await empty_client.fetch_taxonomy_record(taxid)

        assert type(taxonomy) is NCBITaxonomy

        assert {
            path.name for path in empty_client.cache._taxonomy_path.glob("*.json")
        } == snapshot

    @pytest.mark.ncbi
    @pytest.mark.flaky(reruns=3, reruns_delay=3, only_rerun=SERVER_ERRORS)
    @pytest.mark.parametrize("taxid", [1000000000000, 99999999])
    async def test_fetch_taxonomy_from_ncbi_fail(self, taxid: int, empty_client):
        taxonomy = await empty_client.fetch_taxonomy_record(taxid)

        assert taxonomy is None

    @pytest.mark.parametrize("taxid", [438782, 1198450])
    async def test_fetch_taxonomy_from_cache(
        self, taxid: int, scratch_client, snapshot: SnapshotAssertion
    ):
        assert scratch_client.cache.load_taxonomy(taxid)

        taxonomy = await scratch_client.fetch_taxonomy_record(taxid)

        assert type(taxonomy) is NCBITaxonomy

        assert taxonomy == snapshot(name=f"{taxonomy}_validated.json")


@pytest.mark.ncbi
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
@pytest.mark.parametrize(
    "misspelled,expected",
    [
        (
            "Hynchosia yellow mosaic India virus",
            "rhynchosia yellow mosaic india virus",
        ),
        (
            "Angelica bush stunt virus",
            "angelica bushy stunt virus",
        ),
    ],
)
async def test_fetch_spelling(misspelled: str, expected: str):
    taxon_name = await NCBIClient.fetch_spelling(name=misspelled)

    assert taxon_name == expected
