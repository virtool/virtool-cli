from collections import namedtuple
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import orjson
from pydantic import BaseModel

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


class MolType(StrEnum):
    """The in vivo molecule type of a sequence.

    Corresponds to Genbank's moltype field
    """

    DNA = "DNA"
    RNA = "RNA"
    TRNA = "tRNA"
    MRNA = "mRNA"
    CRNA = "cRNA"


class Strandedness(StrEnum):
    """Strandedness of a molecule, either single or double"""

    SINGLE = "single"
    DOUBLE = "double"


class Topology(StrEnum):
    """Topology of a molecule, either linear or circular"""

    LINEAR = "linear"
    CIRCULAR = "circular"


@dataclass
class Molecule:
    """The strandedness, molecule type and topology of this OTU.

    Corresponds to all sequences found in this OTU.
    """

    strandedness: Strandedness
    type: MolType
    topology: Topology


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


IsolateNameKey = namedtuple("IsolateNameKey", ["type", "value"])


class IsolateName(BaseModel):
    """Represents a sub-species categorization name for sequences.

    Can be typed as an isolate, a strain or a clone.
    """

    type: IsolateNameType
    """The type of sub-species categorization"""

    value: str
    """The name of this subcategory"""

    @property
    def frozen(self):
        return IsolateNameKey(type=self.type, value=self.value)


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
