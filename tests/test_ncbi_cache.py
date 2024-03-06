import json

import pytest

from virtool_cli.ncbi.cache import NCBICache


@pytest.fixture()
def empty_cache_path(tmp_path):
    return tmp_path / "test_cache"


def get_test_record(accession: str, cache_example_path):
    with open(cache_example_path / "nuccore" / f"{accession}.json", "r") as f:
        return json.load(f)


def get_test_taxonomy(taxon_id: int, cache_example_path):
    with open(cache_example_path / "taxonomy" / f"{taxon_id}.json", "r") as f:
        return json.load(f)


def test_cache_init(empty_cache_path):
    cache = NCBICache(path=empty_cache_path)

    assert cache.nuccore.is_dir()

    assert cache.taxonomy.is_dir()

    assert cache.nuccore.exists()

    assert cache.taxonomy.exists()


def test_cache_clear(cache_scratch_path):
    cache = NCBICache(path=cache_scratch_path)

    cache.clear()

    assert not list(cache.nuccore.glob("*.json"))

    assert not list(cache.taxonomy.glob("*.json"))


@pytest.mark.parametrize(
    "accessions",
    (
        ["AB017504", "MH200607", "MK431779", "NC_003355"],
        ["NC_036587", "MT240513", "MT240490"],
    ),
)
class TestCacheNuccoreOperations:
    def test_cache_nuccore_load_records(self, accessions, cache_scratch_path):
        scratch_cache = NCBICache(cache_scratch_path)

        records = scratch_cache.load_nuccore_records(accessions)

        print(records)

        assert type(records) is list

        for record in records:
            assert type(record) is dict

    def test_cache_nuccore_cache_records(
        self, accessions, cache_example_path, empty_cache_path
    ):
        fresh_cache = NCBICache(empty_cache_path)

        records = []
        for accession in accessions:
            record = get_test_record(accession, cache_example_path)
            records.append(record)

        fresh_cache.cache_nuccore_records(records)

        for accession in accessions:
            assert (fresh_cache.nuccore / f"{accession}.json").exists()


@pytest.mark.parametrize("taxid", (270478, 438782, 1198450))
class TestCacheTaxonomyOperations:
    def test_cache_taxonomy_load(self, taxid, cache_scratch_path):
        scratch_cache = NCBICache(cache_scratch_path)

        taxonomy = scratch_cache.load_taxonomy(taxid)

        assert type(taxonomy) is dict

    def test_cache_taxonomy_cache(self, taxid, cache_example_path, empty_cache_path):
        taxonomy = get_test_taxonomy(taxid, cache_example_path)

        fresh_cache = NCBICache(empty_cache_path)

        fresh_cache.cache_taxonomy(taxonomy, taxid)

        assert (fresh_cache.taxonomy / f"{taxid}.json").exists()
