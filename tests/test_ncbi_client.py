import pytest
import asyncio
import socket

from urllib.error import HTTPError

from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.ncbi.model import NCBINuccore, NCBISource, NCBITaxonomy


@pytest.fixture()
def blocked_socket():
    blocked_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    blocked_socket.setblocking(True)
    return blocked_socket


@pytest.fixture()
def scratch_client(cache_scratch_path):
    return NCBIClient(cache_scratch_path)


@pytest.fixture()
def empty_client(tmp_path):
    return NCBIClient(NCBICache(tmp_path / "clean_repo").path)


class TestClientFetchGenbank:
    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    async def test_fetch_genbank_records_from_ncbi(self, accessions, scratch_client):
        clean_records = await scratch_client.fetch_genbank_records(
            accessions=accessions, cache_results=True, use_cached=False
        )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.source) is NCBISource

            assert (
                scratch_client.cache.load_nuccore_record(record.accession) is not None
            )

    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    async def test_fetch_genbank_record_from_cache(
        self, accessions, cache_scratch_path, blocked_socket
    ):
        client = NCBIClient(cache_scratch_path)

        with blocked_socket:
            clean_records = await client.fetch_genbank_records(
                accessions=accessions, cache_results=False, use_cached=True
            )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.source) is NCBISource

    @pytest.mark.parametrize(
        "accessions", [["AB017503", "AB017504", "MH200607", "MK431779", "NC_003355"]]
    )
    async def test_fetch_partially_cached_genbank_records(
        self, accessions, cache_scratch_path
    ):
        client = NCBIClient(cache_scratch_path)

        assert next(client.cache.nuccore.glob("*.json"))

        clean_records = await client.fetch_genbank_records(
            accessions=accessions, cache_results=False, use_cached=True
        )

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.source) is NCBISource

    @pytest.mark.asyncio
    async def test_fetch_accessions_fail(self, scratch_client):
        false_accessions = ["friday", "paella", "111"]

        records = await scratch_client.fetch_genbank_records(false_accessions)

        assert not records


@pytest.mark.parametrize(
    "accessions",
    [
        ["AB017504", "MH200607", "MK431779", "NC_003355"],
        ["NC_036587", "MT240513", "MT240490"],
    ],
)
class TestClientFetchRawGenbank:
    @pytest.mark.asyncio
    async def test_fetch_raw_via_accessions(self, accessions):
        records = await NCBIClient.fetch_unvalidated_genbank_records(accessions)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_raw_via_accessions_partial(self, accessions):
        partial_accession_list = accessions
        partial_accession_list[0] = partial_accession_list[0][:3]

        records = await NCBIClient.fetch_unvalidated_genbank_records(accessions)

        assert len(records) == len(accessions) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)


@pytest.mark.asyncio
@pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
async def test_fetch_records_by_taxid(taxid, empty_client):
    with pytest.raises(StopIteration):
        next(empty_client.cache.nuccore.glob("*.json"))

    records = await empty_client.link_from_taxid_and_fetch(taxid, cache_results=True)

    assert records

    for record in records:
        assert type(record) is NCBINuccore

        assert empty_client.cache.load_nuccore_record(record.accession) is not None


class TestClientFetchTaxonomy:
    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
    async def test_fetch_taxonomy_from_ncbi(self, taxid, empty_client):
        for attempt in range(3):
            try:
                taxonomy = await empty_client.fetch_taxonomy(
                    taxid, use_cached=False, cache_results=True
                )

                assert type(taxonomy) is NCBITaxonomy

                assert empty_client.cache.load_taxonomy(taxid) is not None

                return
            except HTTPError:
                await asyncio.sleep(3)

        assert False

    @pytest.mark.asyncio
    @pytest.mark.parametrize("taxid", [438782, 1198450])
    async def test_fetch_taxonomy_from_cache(
        self, taxid, scratch_client, blocked_socket
    ):
        assert scratch_client.cache.load_taxonomy(taxid)

        with blocked_socket:
            taxonomy = await scratch_client.fetch_taxonomy(taxid, use_cached=True)

        assert type(taxonomy) is NCBITaxonomy


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "name,taxid",
    [
        ("Abaca bunchy top virus", 438782),
        ("Rhynchosia golden mosaic virus", 117198),
    ],
)
async def test_fetch_taxonomy_by_name(name, taxid):
    assert await NCBIClient.fetch_taxonomy_id_by_name(name) == taxid


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
async def test_check_spelling(misspelled, expected):
    taxon_name = await NCBIClient.check_spelling(name=misspelled)

    assert taxon_name.lower() == expected.lower()
