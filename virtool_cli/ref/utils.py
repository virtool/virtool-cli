from dataclasses import dataclass
from enum import StrEnum
import re
from pathlib import Path
from uuid import UUID

import orjson

from virtool_cli.legacy.models import LegacyIsolateSource, LegacySourceType
from virtool_cli.ncbi.model import NCBIGenbank


class DataType(StrEnum):
    BARCODE = "barcode"
    GENOME = "genome"


def format_json(path: Path):
    """Format the JSON file at `path` in place."""
    with path.open("rb") as f:
        data = orjson.loads(f.read())

    with path.open("wb") as f:
        f.write(orjson.dumps(data, option=orjson.OPT_INDENT_2))


def pad_zeroes(number: int) -> str:
    """Pad a number with zeroes to make it 8 characters long.

    :param number: the number to pad
    :return: the padded number
    """
    if number > 99999999:
        raise ValueError("Number is too large to pad")

    return str(number).zfill(8)


class IsolateNameType(StrEnum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    REFSEQ = "refseq"

    def __repr__(self):
        return self.value


isolate_name_pattern = re.compile(r"IsolateName\(type=(\w+), value=[\'\"](.+)[\'\"]\)")


@dataclass(frozen=True)
class IsolateName:
    """Represents a sub-species categorization name for sequences.

    Can be typed as an isolate, a strain or a clone.
    """

    type: IsolateNameType
    """The type of sub-species categorization"""

    value: str
    """The name of this subcategory"""

    @classmethod
    def from_dict(cls, raw: dict):
        return IsolateName(type=IsolateNameType(raw["type"]), value=raw["value"])

    @classmethod
    def from_string(cls, raw: str):
        parts = re.match(isolate_name_pattern, raw).groups()

        return IsolateName(type=IsolateNameType(parts[0]), value=parts[1])


def extract_isolate_source(
    genbank_records: list[NCBIGenbank],
) -> LegacyIsolateSource:
    """Extract a legacy isolate source from a set of Genbank records associated with the
    isolate.
    """
    for record in genbank_records:
        if record.source.isolate:
            return LegacyIsolateSource(
                name=record.source.isolate,
                type=LegacySourceType.ISOLATE,
            )

        if record.source.strain:
            return LegacyIsolateSource(
                name=record.source.strain,
                type=LegacySourceType.STRAIN,
            )

        if record.source.clone:
            return LegacyIsolateSource(
                name=record.source.clone,
                type=LegacySourceType.CLONE,
            )

    accessions = sorted(
        (record.accession for record in genbank_records if record.accession),
        key=lambda x: int(x.replace("NC_", "").replace(".1", "")),
    )

    return LegacyIsolateSource(
        name=accessions[0].upper(),
        type=LegacySourceType.GENBANK,
    )


def parse_uuid(raw: str) -> UUID | None:
    """Returns UUID or None without erroring"""
    try:
        return UUID(raw)
    except ValueError:
        return None
