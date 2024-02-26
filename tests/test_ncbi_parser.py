import json
from pathlib import Path

import pytest

from virtool_cli.ncbi.parser import RecordParser, TaxonomyParser
from virtool_cli.ncbi.error import NCBIParseError
from virtool_cli.ncbi.model import NCBIAccession, NCBISource, NCBITaxonomy, NCBIRank


@pytest.fixture
def test_records_path(test_files_path: Path):
    return test_files_path / "cache_test" / "nuccore"


@pytest.fixture
def test_taxonomy_path(test_files_path: Path):
    return test_files_path / "cache_test" / "taxonomy"


def get_test_record_set(filename: str, test_records_path):
    with open(test_records_path / f"{filename}.json", "r") as f:
        return json.load(f)


def get_test_taxonomy(taxon_id: int, test_taxonomy_path):
    with open(test_taxonomy_path / f"{taxon_id}.json", "r") as f:
        return json.load(f)


class TestRecordParser:
    @pytest.mark.parametrize("record_otu", ["0bfdb8bc", "2ksa35mn"])
    def test_init_parser_success(self, record_otu, test_records_path):
        record_set = get_test_record_set(record_otu, test_records_path)

        for record in record_set:
            record_parser = RecordParser(record)
            assert type(record_parser) is RecordParser

    @pytest.mark.parametrize(
        "bad_data",
        [
            {},
            {"id": "AJ2534F.1"},
            {"GBSeq_accession-version": "AJ2534F.1"},
        ],
    )
    def test_init_parser_fail(self, bad_data):
        with pytest.raises(NCBIParseError):
            RecordParser(data=bad_data)

        try:
            RecordParser(data=bad_data)
        except NCBIParseError as e:
            key_complete = (
                "GBSeq_accession-version",
                "GBSeq_sequence",
                "GBSeq_definition",
                "GBSeq_feature-table",
            ) in e.keys
            assert key_complete is False

    @pytest.mark.parametrize("record_otu", ["0bfdb8bc", "2ksa35mn"])
    def test_read_success(self, record_otu, test_records_path):
        record_set = get_test_record_set(record_otu, test_records_path)

        for record in record_set:
            record_parser = RecordParser(record)

            assert type(record_parser) is RecordParser

            primary, source = record_parser.read()

            assert type(primary) is NCBIAccession

            assert type(source) is NCBISource

            assert type(primary.accession) is str

            assert type(source.taxid) is int


class TestTaxonomyParser:
    @pytest.mark.parametrize("taxon_id", [438782, 1198450, 270478])
    def test_parse_taxonomy(self, taxon_id, test_taxonomy_path):
        taxonomy = get_test_taxonomy(taxon_id, test_taxonomy_path)

        taxonomy_parser = TaxonomyParser(taxonomy)

        assert type(taxonomy_parser) is TaxonomyParser

        assert type(taxonomy_parser.get_taxon_id()) is int

        assert type(taxonomy_parser.get_lineage()) is list

    @pytest.mark.parametrize("taxon_id", [438782, 1198450, 270478])
    def test_read_taxonomy(self, taxon_id, test_taxonomy_path):
        raw_taxonomy = get_test_taxonomy(taxon_id, test_taxonomy_path)

        taxonomy_parser = TaxonomyParser(raw_taxonomy)

        parsed_taxonomy = taxonomy_parser.read()

        assert type(parsed_taxonomy) is NCBITaxonomy

        assert type(parsed_taxonomy.id) is int

        assert type(parsed_taxonomy.rank) is NCBIRank

        assert type(parsed_taxonomy.lineage) is list

        assert parsed_taxonomy.accessions == []
