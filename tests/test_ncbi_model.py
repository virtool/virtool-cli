import pytest

from pydantic import ValidationError
from syrupy import SnapshotAssertion

from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.ncbi.model import NCBIGenbank, NCBISource, NCBILineage, NCBITaxonomy


@pytest.fixture()
def scratch_cache(cache_scratch_path):
    return NCBICache(cache_scratch_path)


@pytest.mark.parametrize(
    "accession", ["AB017504", "MH200607", "NC_036587", "MT240513", "NC_015504"]
)
class TestParseGenbank:
    def test_parse_genbank_source(
        self, accession, scratch_cache, snapshot: SnapshotAssertion
    ):
        record = scratch_cache.load_genbank_record(accession)

        source = NCBIGenbank.create_source(record["GBSeq_feature-table"])

        assert source.model_dump() == snapshot

    def test_parse_genbank_record(
        self, accession, scratch_cache, snapshot: SnapshotAssertion
    ):
        record = scratch_cache.load_genbank_record(accession)

        validated_record = NCBIGenbank(**record)

        assert validated_record.model_dump() == snapshot

    def test_parse_genbank_source_taxid(
        self, accession, scratch_cache, snapshot: SnapshotAssertion
    ):
        record = scratch_cache.load_genbank_record(accession)

        db_xref = None
        for feature in record["GBSeq_feature-table"]:
            if feature["GBFeature_key"] == "source":
                for qual in feature["GBFeature_quals"]:
                    if qual["GBQualifier_name"] == "db_xref":
                        db_xref = qual["GBQualifier_value"]
                        break

        assert NCBISource.db_xref_to_taxid(db_xref) == snapshot


class TestTaxonomyParse:
    @pytest.mark.parametrize("taxid", [270478, 438782, 1198450])
    def test_parse_taxonomy_record_single(
        self, taxid: int, scratch_cache, snapshot: SnapshotAssertion
    ):
        record = scratch_cache.load_taxonomy(taxid)

        validated_taxonomy = NCBITaxonomy(**record)

        assert type(validated_taxonomy) is NCBITaxonomy

        assert validated_taxonomy.model_dump() == snapshot

    @pytest.mark.parametrize("taxid, rank", [(1016856, "isolate")])
    def test_parse_taxonomy_record_rank_addendum(
        self, taxid: int, rank: str, scratch_cache, snapshot: SnapshotAssertion
    ):
        record = scratch_cache.load_taxonomy(taxid)

        with pytest.raises(ValidationError):
            assert NCBITaxonomy(**record)

        validated_taxonomy = NCBITaxonomy(rank=rank, **record)

        assert validated_taxonomy.model_dump() == snapshot


@pytest.mark.parametrize(
    "lineage_data",
    [
        {"TaxId": "2732397", "ScientificName": "Pararnavirae", "Rank": "kingdom"},
        {"TaxId": "2732409", "ScientificName": "Artverviricota", "Rank": "phylum"},
    ],
)
def test_create_lineage_item_alias(lineage_data, snapshot: SnapshotAssertion):
    lineage_data_clean = {
        "id": lineage_data["TaxId"],
        "name": lineage_data["ScientificName"],
        "rank": lineage_data["Rank"],
    }

    lineage_item_native = NCBILineage(**lineage_data)
    lineage_item_clean = NCBILineage(**lineage_data_clean)

    assert lineage_item_native == lineage_item_clean

    assert lineage_item_clean == snapshot
