from enum import StrEnum

from virtool_cli.ncbi.error import NCBIParseError
from virtool_cli.ncbi.model import (
    NCBIAccession,
    NCBISourceType,
    NCBISource,
    NCBITaxonomy,
    NCBIRank,
)


class GBSeq(StrEnum):
    ACCESSION = "GBSeq_accession-version"
    DEFINITION = "GBSeq_definition"
    SEQUENCE = "GBSeq_sequence"
    LENGTH = "GBSeq_length"
    COMMENT = "GBSeq_comment"
    FEATURE_TABLE = "GBSeq_feature-table"


class RecordParser:
    def __init__(self, data: dict):
        self.raw = data

        try:
            self.source_parser = SourceParser(self._get_source_table())
        except KeyError as e:
            raise NCBIParseError(keys=data, message=str(e))

    def read(self) -> tuple[NCBIAccession, NCBISource]:
        try:
            return self.read_primary(), self.read_source()
        except KeyError as e:
            raise NCBIParseError(keys=self.raw.keys, message=str(e))

    def read_primary(self) -> NCBIAccession:
        return NCBIAccession(
            accession=self.raw[GBSeq.ACCESSION],
            definition=self.raw[GBSeq.DEFINITION],
            sequence=self.raw[GBSeq.SEQUENCE],
            comment=self.raw.get(GBSeq.COMMENT, ""),
        )

    def read_source(self) -> NCBISource:
        return self.source_parser.read()

    @staticmethod
    def get_accession(record: dict) -> str:
        return record[GBSeq.ACCESSION]

    @staticmethod
    def get_definition(record: dict) -> str:
        return record[GBSeq.DEFINITION]

    @staticmethod
    def get_length(record: dict) -> int:
        return record[GBSeq.LENGTH]

    @staticmethod
    def get_comment(record: dict) -> str | None:
        return record.get(GBSeq.COMMENT, "")

    @staticmethod
    def get_sequence(record: dict) -> str:
        return record[GBSeq.SEQUENCE]

    @staticmethod
    def get_feature_table(record: dict) -> list[dict]:
        return record[GBSeq.FEATURE_TABLE]

    def get_taxid(self) -> int:
        return self.source_parser.get_taxid()

    def get_host(self) -> str:
        return self.source_parser.get_host()

    def get_segment(self) -> str:
        return self.source_parser.get_segment()

    def _get_source_table(self) -> dict | None:
        for feature in self.raw[GBSeq.FEATURE_TABLE]:
            if feature["GBFeature_key"] == "source":
                return feature

        raise NCBIParseError(
            keys=self._get_feature_table_keys(),
            message="Feature table does not contain source data",
        )

    def _get_feature_table_keys(self):
        keys = []
        for feature in self.raw[GBSeq.FEATURE_TABLE]:
            keys.append(feature["GBFeature_key"])

        return keys


class FeatureParser:
    def __init__(self, data: dict):
        self.raw = data

        self.dictionary = self.to_dict(data)

    def read(self):
        """
        Placeholder. One-step parse into an object.
        """

    @staticmethod
    def to_dict(feature: dict) -> dict:
        """Converts the feature table format to a dict"""
        print(feature)

        qualifier_dict = {}
        for qualifier in feature["GBFeature_quals"]:
            qual_name = qualifier["GBQualifier_name"]
            qual_value = qualifier["GBQualifier_value"]
            qualifier_dict[qual_name] = qual_value

        return qualifier_dict

    def _get_error(self, message: str = ""):
        return NCBIParseError(keys=self.dictionary.keys, message=message)


class SourceParser(FeatureParser):
    """Parses the source feature in a GenBank record"""

    def __init__(self, data: dict):
        super().__init__(data)

    def read(self) -> NCBISource:
        try:
            source_type = self.get_source_type()
        except KeyError:
            raise NCBIParseError(
                keys=self.dictionary.keys,
                message="Not enough data in this source table",
            )

        source_name = self._get_source_name(source_type)

        return NCBISource(
            type=source_type,
            name=source_name,
            host=self.get_host(),
            segment=self.get_segment(),
            taxid=self.get_taxid(),
        )

    def get_host(self) -> str:
        return self.dictionary.get("host", "")

    def get_segment(self) -> str:
        return self.dictionary.get("segment", "")

    def get_taxid(self) -> int | None:
        value = self.dictionary["db_xref"]

        key, taxid = value.split(":")

        if key == "taxon":
            return taxid

        return None

    def get_source_type(self) -> NCBISourceType:
        for key in self.dictionary:
            try:
                return NCBISourceType(key)
            except ValueError:
                continue

        raise NCBIParseError(
            keys=self.dictionary.keys, message="Missing source type qualifier"
        )

    def get_source_name(self):
        source_type = self.get_source_type()

        return self._get_source_name(source_type)

    def _get_source_name(self, source_type) -> str | None:
        try:
            return self.dictionary[source_type.value]
        except KeyError:
            raise NCBIParseError(
                keys=self.dictionary.keys,
                message="Not enough data in this source table",
            )


class TaxonomyParser:
    def __init__(self, data):
        self.raw = data

    def read(self):
        return NCBITaxonomy(
            id=int(self.raw["TaxId"]),
            accessions=[],
            lineage=self.get_lineage(),
            rank=self.get_rank(),
        )

    def get_taxon_id(self) -> int:
        return int(self.raw["TaxId"])

    def get_name(self) -> str:
        return self.raw["ScientificName"]

    def get_lineage(self) -> list[dict] | None:
        lineage = []
        for level in self.raw["LineageEx"]:
            lineage.append(level["ScientificName"])

        return lineage

    def get_rank(self) -> NCBIRank | None:
        try:
            return NCBIRank(self.raw["Rank"])
        except Exception as e:
            raise e


class CDSParser(FeatureParser):
    def __init__(self, data):
        super().__init__(data)

    def get_codon_start(self):
        return self.dictionary["codon_start"]

    def get_transl_table(self):
        return self.dictionary["transl_table"]

    def get_product(self):
        return self.dictionary["product"]

    def get_protein_id(self):
        return self.dictionary["protein_id"]

    def get_translation(self):
        return self.dictionary["translation"]
