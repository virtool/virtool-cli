from typing import Annotated, Dict

from pydantic import (
    BaseModel,
    Field,
    TypeAdapter,
    UUID4,
    field_validator,
)

from virtool_cli.utils.models import Molecule
from virtool_cli.ref.utils import IsolateName, IsolateNameType


class OTUSnapshotSequence(BaseModel):
    """Represents a sequence in a Virtool reference repository."""

    id: UUID4
    """The sequence id."""

    accession: str
    """The sequence accession."""

    definition: str
    """The sequence definition."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the sequence was not migrated from a legacy repository, this will be `None`.
    """

    sequence: str
    """The sequence."""

    segment: str
    """The sequence segment."""


class OTUSnapshotIsolate(BaseModel):
    id: UUID4
    """The isolate ID."""

    name: IsolateName
    """The isolate's source name metadata."""

    legacy_id: str | None = None
    """A string based ID carried over from a legacy Virtool reference repository."""

    @field_validator("name", mode="before")
    @classmethod
    def convert_isolate_name(cls, raw: dict) -> IsolateName:
        """Takes a dictionary and converts to IsolateName."""
        return IsolateName(type=IsolateNameType(raw["type"]), value=raw["value"])


class OTUSnapshotOTU(BaseModel):
    id: UUID4
    """The OTU id."""

    taxid: int
    """The NCBI Taxonomy id for this OTU."""

    name: str
    """The name of the OTU (eg. TMV for Tobacco mosaic virus)"""

    acronym: str = ""
    """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository."""

    molecule: Annotated[Molecule | None, Field()] = None
    """The molecule of this OTU"""

    segment_schema: Annotated[list | None, Field(alias="schema")] = None

    repr_isolate: UUID4 | None = None


class OTUSnapshotToCIsolate(BaseModel):
    id: UUID4
    accessions: dict[str, UUID4]


toc_adapter = TypeAdapter(Dict[str, OTUSnapshotToCIsolate])
