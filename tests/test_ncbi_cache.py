import json

import pytest

from virtool_cli.ncbi.cache import NCBICache


@pytest.fixture()
def empty_cache_path(tmp_path):
    return tmp_path / "test_cache"


def get_test_record(accession: str, cache_example_path) -> dict:
    with open(cache_example_path / "nuccore" / f"{accession}.json", "r") as f:
        return json.load(f)


def get_test_taxonomy(taxon_id: int, cache_example_path) -> dict:
    with open(cache_example_path / "taxonomy" / f"{taxon_id}.json", "r") as f:
        return json.load(f)


def test_cache_init(empty_cache_path):
    assert not empty_cache_path.exists()

    cache = NCBICache(path=empty_cache_path)

    assert cache.nuccore.is_dir()

    assert cache.taxonomy.is_dir()

    assert cache.nuccore.exists()

    assert cache.taxonomy.exists()


def test_cache_clear(cache_scratch_path):
    cache = NCBICache(path=cache_scratch_path)

    try:
        assert next(cache.nuccore.glob("*.json"))
        assert next(cache.taxonomy.glob("*.json"))
    except StopIteration as exc:
        pytest.fail(f"Cache scratch path was not set up correctly: {exc}")

    cache.clear()

    with pytest.raises(StopIteration):
        assert not next(cache.nuccore.glob("*.json"))

    with pytest.raises(StopIteration):
        assert not next(cache.taxonomy.glob("*.json"))


@pytest.mark.parametrize(
    "accessions",
    (
        ["AB017504", "MH200607", "MK431779", "NC_003355"],
        ["NC_036587", "MT240513", "MT240490"],
    ),
)
class TestCacheNuccoreOperations:
    def test_cache_nuccore_load_record_batch(self, accessions, cache_scratch_path):
        scratch_cache = NCBICache(cache_scratch_path)

        for accession in accessions:
            record = scratch_cache.load_nuccore_record(accession)

            assert type(record) is dict

    def test_cache_nuccore_cache_records(
        self, accessions, cache_example_path, empty_cache_path
    ):
        assert not empty_cache_path.exists()

        cache = NCBICache(empty_cache_path)

        for accession in accessions:
            record = get_test_record(accession, cache_example_path)

            cache.cache_nuccore_record(data=record, accession=accession)

            assert (cache.nuccore / f"{accession}.json").exists()


@pytest.mark.parametrize("fake_accession", ["afjshd", "23222", "wheelhouse"])
def test_cache_nuccore_load_fail(fake_accession, cache_scratch_path):
    scratch_cache = NCBICache(cache_scratch_path)

    assert scratch_cache.load_nuccore_record(fake_accession) is None


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
