import pytest
from pathlib import Path

from virtool_cli.repo.cls import Repo
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.model import NCBINuccore, NCBISource


ACCESSION_LISTS = [
    ["KT390494", "KT390496", "KT390501"],
    ["KY702580", "MZ148028.1", "KF915809"],
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
def accession_list():
    return ["KT390494", "KT390496", "KT390501"]


@pytest.fixture()
def scratch_repo(scratch_path):
    return Repo(scratch_path)


@pytest.mark.parametrize("accession_list", ACCESSION_LISTS)
class TestClientProcureAccessions:
    async def test_procure_accessions(self, accession_list, test_client):
        clean_records = await test_client.procure_accessions(requested=accession_list)

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.accession) is str

            assert type(record.source) is NCBISource

            assert type(record.source.taxid) is int


@pytest.mark.parametrize("taxon_id", [1016856, 429130])
class TestClientProcureFromTaxid:
    @pytest.mark.asyncio
    async def test_procure_from_taxid(self, taxon_id, test_client):
        clean_records = await test_client.procure_from_taxid(taxon_id, use_cached=False)

        assert type(clean_records) is list

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.accession) is str

            assert type(record.source) is NCBISource

            assert type(record.source.taxid) is int

    @pytest.mark.asyncio
    async def test_cache_from_taxid(self, taxon_id, test_client):
        await test_client.cache_from_taxid(taxon_id)

        assert test_client.cache.load_nuccore(taxon_id)


class TestClientProcureUpdates:
    @pytest.mark.asyncio
    @pytest.mark.parametrize(
        "otu_id, use_cached",
        (
            [
                ("0bfdb8bc", True),
                ("0bfdb8bc", False),
                ("4c9ddfb7", True),
                ("4c9ddfb7", False),
                ("d226290f", True),
            ]
        ),
    )
    async def test_procure_updates(self, otu_id, use_cached, scratch_repo):
        otu = scratch_repo.get_otu_by_id(otu_id)
        client = NCBIClient.for_repo(scratch_repo.path)

        clean_records = await client.procure_updates(
            otu_id=otu.id,
            taxid=otu.taxid,
            blocked_accessions=otu.blocked_accessions,
            use_cached=use_cached,
        )

        assert type(clean_records) is list

        for record in clean_records:
            assert type(record) is NCBINuccore

            assert type(record.accession) is str

            assert type(record.source) is NCBISource

            assert type(record.source.taxid) is int

    @pytest.mark.parametrize("otu_id", ["0bfdb8bc", "4c9ddfb7", "d226290f"])
    @pytest.mark.asyncio
    async def test_cache_updates(self, otu_id, scratch_repo):
        otu = scratch_repo.get_otu_by_id(otu_id)

        client = NCBIClient.for_repo(scratch_repo.path)

        await client.cache_updates(
            otu_id=otu.id,
            taxid=otu.taxid,
            blocked_accessions=otu.blocked_accessions,
        )

        assert client.cache.load_nuccore(otu_id)


@pytest.mark.parametrize("accession_list", ACCESSION_LISTS)
class TestClientFetchAccessionSets:
    @pytest.mark.asyncio
    async def test_fetch_accessions(self, accession_list):
        records = await NCBIClient.fetch_by_accessions(accession_list)

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.asyncio
    async def test_fetch_accessions_partial(self, accession_list):
        partial_accession_list = accession_list
        partial_accession_list[0] = partial_accession_list[0][:3]

        records = await NCBIClient.fetch_by_accessions(accession_list)

        assert len(records) == len(accession_list) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)


@pytest.mark.parametrize(
    "requested, blocked",
    [
        (ACCESSION_LISTS[0], None),
        (ACCESSION_LISTS[0], []),
        (ACCESSION_LISTS[0], ACCESSION_LISTS[0][:1]),
        (ACCESSION_LISTS[1], None),
        (ACCESSION_LISTS[1], []),
        (ACCESSION_LISTS[1], ACCESSION_LISTS[0][:2]),
    ],
)
def test_filter_accessions(requested, blocked):
    filtered = NCBIClient.filter_accessions(requested, blocked)

    assert filtered

    if blocked is None:
        return

    for blocked_accession in blocked:
        assert blocked_accession not in filtered


@pytest.mark.asyncio
async def test_fetch_accessions_fail():
    accession_list = ["friday", "paella", "111"]

    assert not await NCBIClient.fetch_by_accessions(accession_list)


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
