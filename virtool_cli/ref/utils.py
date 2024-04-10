from enum import StrEnum
from dataclasses import dataclass


class DataType(StrEnum):
    """The possible reference data types."""

    BARCODE = "barcode"
    GENOME = "genome"


class MolType(StrEnum):
    """The in vivo molecule type of a sequence, corresponds to Genbank's moltype field"""

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

    Corresponds to all sequences found in this OTU."""

    strandedness: Strandedness
    type: MolType
    topology: Topology


def pad_zeroes(number: int) -> str:
    """Pad a number with zeroes to make it 8 characters long.

    :param number: the number to pad
    :return: the padded number
    """
    return str(number).zfill(8)
