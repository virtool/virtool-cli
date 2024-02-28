from dataclasses import dataclass
from enum import StrEnum

from virtool_cli.ncbi.error import NCBIParseError
from virtool_cli.ncbi.model import NCBIAccession, NCBISource, NCBISourceType


class GBSeq(StrEnum):
    ACCESSION = "GBSeq_accession-version"
    DEFINITION = "GBSeq_definition"
    SEQUENCE = "GBSeq_sequence"
    LENGTH = "GBSeq_length"
    COMMENT = "GBSeq_comment"
    FEATURE_TABLE = "GBSeq_feature-table"


@dataclass
class NuccorePacket:
    sequence: NCBIAccession
    source: NCBISource


def parse_nuccore(raw: dict) -> NuccorePacket:
    sequence = NCBIAccession(
        accession=raw[GBSeq.ACCESSION],
        definition=raw[GBSeq.DEFINITION],
        sequence=raw[GBSeq.SEQUENCE],
        comment=raw.get(GBSeq.COMMENT, ""),
    )

    source = parse_source(feature_table_to_dict(get_source_table(raw)))

    return NuccorePacket(sequence, source)


def parse_source(source_dict: dict) -> NCBISource:
    try:
        source_type = get_source_type(source_dict)
        source_name = source_dict[source_type]
    except KeyError:
        raise NCBIParseError(
            keys=source_dict.keys,
            message="Not enough data in this source table",
        )

    return NCBISource(
        type=source_type,
        name=source_name,
        host=source_dict.get("host", ""),
        segment=source_dict.get("segment", ""),
        taxid=parse_taxid(source_dict),
    )


def get_source_type(source_qualifiers: dict) -> NCBISourceType:
    """Return a NCBISourceType, prioritizing ISOLATE over other options."""
    if NCBISourceType.ISOLATE in source_qualifiers:
        return NCBISourceType.ISOLATE

    if NCBISourceType.STRAIN in source_qualifiers:
        return NCBISourceType.STRAIN

    if NCBISourceType.CLONE in source_qualifiers:
        return NCBISourceType.CLONE

    if NCBISourceType.GENOTYPE in source_qualifiers:
        return NCBISourceType.GENOTYPE

    raise NCBIParseError(
        keys=list(source_qualifiers.keys()), message="Missing source type qualifier"
    )


def get_source_table(raw) -> dict | None:
    for feature in raw[GBSeq.FEATURE_TABLE]:
        if feature["GBFeature_key"] == "source":
            return feature

    raise NCBIParseError(
        keys=get_feature_table_keys(raw),
        message="Feature table does not contain source data",
    )


def get_feature_table_keys(raw):
    keys = []
    for feature in raw[GBSeq.FEATURE_TABLE]:
        keys.append(feature["GBFeature_key"])

    return keys


def feature_table_to_dict(feature: dict) -> dict:
    """Converts the feature table format to a dict"""
    qualifier_dict = {}
    for qualifier in feature["GBFeature_quals"]:
        qual_name = qualifier["GBQualifier_name"]
        qual_value = qualifier["GBQualifier_value"]
        qualifier_dict[qual_name] = qual_value

    return qualifier_dict


def parse_taxid(source_feature) -> int | None:
    value = source_feature["db_xref"]

    key, taxid = value.split(":")

    return taxid
