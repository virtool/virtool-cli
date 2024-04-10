from enum import StrEnum

from pydantic import BaseModel, field_validator


class SubSpeciesBin(StrEnum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    VARIANT = "variant"


class IsolateName(BaseModel):
    """Represents a sub-species subcategorization label for sequences"""

    type: SubSpeciesBin
    """The nature of of the subcategory."""

    value: str
    """The name of the subcategory"""

    @field_validator("type", mode="before")
    @classmethod
    def validate_type(cls, raw: str):
        return SubSpeciesBin(raw)
