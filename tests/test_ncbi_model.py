import pytest
from syrupy import SnapshotAssertion

from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.ncbi.model import NCBINuccore, NCBISource, NCBILineage


@pytest.fixture()
def scratch_cache(cache_scratch_path):
    return NCBICache(cache_scratch_path)


@pytest.mark.parametrize("accession", ["AB017504", "MH200607", "NC_036587", "MT240513"])
class TestAccessionParse:
    def test_parse_source(self, accession, scratch_cache, snapshot: SnapshotAssertion):
        record = scratch_cache.load_nuccore_record(accession)

        source = NCBINuccore.create_source(record["GBSeq_feature-table"])

        assert source == snapshot

    def test_parse_nuccore(self, accession, scratch_cache, snapshot: SnapshotAssertion):
        record = scratch_cache.load_nuccore_record(accession)

        validated_record = NCBINuccore(**record)

        assert validated_record == snapshot

    def test_parse_source_taxid(
        self, accession, scratch_cache, snapshot: SnapshotAssertion
    ):
        record = scratch_cache.load_nuccore_record(accession)

        db_xref = None
        for feature in record["GBSeq_feature-table"]:
            if feature["GBFeature_key"] == "source":
                for qual in feature["GBFeature_quals"]:
                    if qual["GBQualifier_name"] == "db_xref":
                        db_xref = qual["GBQualifier_value"]
                        break

        assert NCBISource.db_xref_to_taxid(db_xref) == snapshot


@pytest.mark.parametrize(
    "lineage_data",
    [
        {"TaxId": "2732397", "ScientificName": "Pararnavirae", "Rank": "kingdom"},
        {"TaxId": "2732409", "ScientificName": "Artverviricota", "Rank": "phylum"},
    ],
)
def test_create_lineage_item_alias(lineage_data, snapshot: SnapshotAssertion):
    lineage_item_native = NCBILineage(**lineage_data)
    lineage_item_clean = NCBILineage(
        id=lineage_data["TaxId"],
        name=lineage_data["ScientificName"],
        rank=lineage_data["Rank"],
    )

    assert lineage_item_native == lineage_item_clean

    assert lineage_item_clean == snapshot
