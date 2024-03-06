import pytest
import json
from pathlib import Path

from virtool_cli.repo.cls import Repo
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.model import NCBINuccore, NCBISource, NCBITaxonomy


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
    @staticmethod
    def check_record(record, cache):
        assert type(record) is NCBINuccore

        assert type(record.accession) is str

        assert type(record.source) is NCBISource

        assert type(record.source.taxid) is int

        assert cache.load_nuccore_record(record.accession.split(".")[0])

    @pytest.mark.parametrize("accessions", ACCESSION_LISTS)
    async def test_fetch_accessions_from_ncbi(self, accessions, test_client):
        clean_records = await test_client.fetch_accessions(
            accessions=accessions, cache_results=True, use_cached=False
        )

        assert clean_records

        for record in clean_records:
            self.check_record(record, test_client.cache)

    @pytest.mark.parametrize("accessions", ACCESSION_LISTS)
    async def test_fetch_accessions_from_cache(self, accessions, test_client):
        clean_records = await test_client.fetch_accessions(
            accessions=accessions, cache_results=True, use_cached=True
        )

        assert clean_records

        for record in clean_records:
            self.check_record(record, test_client.cache)

    @pytest.mark.asyncio
    async def test_fetch_accessions_fail(self, test_client):
        false_accessions = ["friday", "paella", "111"]

        records = await test_client.fetch_accessions(false_accessions)

        assert not records


@pytest.mark.parametrize("accessions", ACCESSION_LISTS)
class TestClientFetchRawAccessions:
    @pytest.mark.asyncio
    async def test_fetch_raw_via_accessions(self, accessions):
        records = await NCBIClient.fetch_unvalidated_accessions(accessions)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_raw_via_accessions_partial(self, accessions):
        partial_accession_list = accessions
        partial_accession_list[0] = partial_accession_list[0][:3]

        records = await NCBIClient.fetch_unvalidated_accessions(accessions)

        assert len(records) == len(accessions) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)


@pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
class TestClientFetchTaxonomy:
    @pytest.mark.asyncio
    async def test_fetch_taxonomy_from_cache(self, taxid, test_client):
        taxonomy = await test_client.fetch_taxonomy(taxid, use_cached=True)

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
