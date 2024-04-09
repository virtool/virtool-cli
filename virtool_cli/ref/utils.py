from enum import StrEnum
from dataclasses import dataclass


class DataType(StrEnum):
    """The possible reference data types."""

    BARCODE = "barcode"
    GENOME = "genome"


class MolType(StrEnum):
    DNA = "DNA"
    RNA = "RNA"
    TRNA = "tRNA"
    MRNA = "mRNA"
    CRNA = "cRNA"


class Strandedness(StrEnum):
    SINGLE = "single"
    DOUBLE = "double"


class Topology(StrEnum):
    LINEAR = "linear"
    CIRCULAR = "circular"


@dataclass
class Molecule:
    strandedness: Strandedness
    type: MolType
    topology: Topology


def pad_zeroes(number: int) -> str:
    """Pad a number with zeroes to make it 8 characters long.

    :param number: the number to pad
    :return: the padded number
    """
    return str(number).zfill(8)
