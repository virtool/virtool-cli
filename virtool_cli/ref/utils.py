from enum import StrEnum


class DataType(StrEnum):
    """The possible reference data types."""

    BARCODE = "barcode"
    GENOME = "genome"


def pad_zeroes(number: int) -> str:
    """Pad a number with zeroes to make it 8 characters long.

    :param number: the number to pad
    :return: the padded number
    """
    return str(number).zfill(8)
