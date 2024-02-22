import json

import pytest

from virtool_cli.ncbi.cache import NCBICache


@pytest.fixture()
def empty_cache_path(tmp_path):
    return tmp_path / "test_cache"


def get_test_record_set(filename: str, cache_example_path):
    with open(cache_example_path / "nuccore" / f"{filename}.json", "r") as f:
        return json.load(f)


def get_test_taxonomy(taxon_id: int, cache_example_path):
    with open(cache_example_path / "taxonomy" / f"{taxon_id}.json", "r") as f:
        return json.load(f)


class TestCache:
    def test_cache_init(self, empty_cache_path):
        cache = NCBICache(path=empty_cache_path)

        assert cache.nuccore.is_dir()

        assert cache.taxonomy.is_dir()

        assert cache.nuccore.exists()

        assert cache.taxonomy.exists()

    def test_cache_clear(self, cache_scratch_path):
        cache = NCBICache(path=cache_scratch_path)

        cache.clear()

        assert not list(cache.nuccore.glob("*.json"))

        assert not list(cache.taxonomy.glob("*.json"))

    @pytest.mark.parametrize("record_otu", ["0bfdb8bc", "2ksa35mn"])
    def test_cache_nuccore_load(self, record_otu, cache_scratch_path):
        scratch_cache = NCBICache(cache_scratch_path)

        records = scratch_cache.load_records(record_otu)

        assert type(records) is list

        for record in records:
            assert type(record) is dict

    @pytest.mark.parametrize("record_otu", ["0bfdb8bc", "2ksa35mn"])
    def test_cache_nuccore_cache(
        self, record_otu, cache_example_path, empty_cache_path
    ):
        fresh_cache = NCBICache(empty_cache_path)

        records = get_test_record_set(record_otu, cache_example_path)

        fresh_cache.cache_records(records, filestem=record_otu)

        assert (fresh_cache.nuccore / f"{record_otu}.json").exists()

    @pytest.mark.parametrize("record_otu", ["0bfdb8bc", "2ksa35mn"])
    def test_cache_taxonomy_load(self, record_otu, cache_scratch_path):
        scratch_cache = NCBICache(cache_scratch_path)

        records = scratch_cache.load_records(record_otu)

        assert type(records) is list

        for record in records:
            assert type(record) is dict

    @pytest.mark.parametrize("record_otu", ["0bfdb8bc", "2ksa35mn"])
    def test_cache_taxonomy_cache(
        self, record_otu, cache_example_path, empty_cache_path
    ):
        fresh_cache = NCBICache(empty_cache_path)

        records = get_test_record_set(record_otu, cache_example_path)

        fresh_cache.cache_records(records, filestem=record_otu)

        assert (fresh_cache.nuccore / f"{record_otu}.json").exists()
