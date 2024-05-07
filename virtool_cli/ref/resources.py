import datetime
from dataclasses import dataclass
from pydantic import BaseModel, field_validator
from uuid import UUID

from virtool_cli.ref.utils import DataType, IsolateName, Molecule


class RepoMeta(BaseModel):
    """Represents the metadata for a Virtool reference repository."""

    id: UUID
    """The repository id."""

    created_at: datetime.datetime
    """The date and time the repository was created."""

    data_type: DataType
    """The repository data type."""

    name: str
    """The repository name."""

    organism: str
    """The repository organism."""

    @field_validator("data_type", mode="before")
    @classmethod
    def datatype_conversion(cls, raw: str) -> DataType:
        return DataType(raw)


@dataclass
class EventSourcedRepoSequence:
    """Represents a sequence in a Virtool reference repository."""

    id: UUID
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

    def dict(self):
        return {
            "id": self.id,
            "accession": self.accession,
            "definition": self.definition,
            "legacy_id": self.legacy_id,
            "sequence": self.sequence,
            "segment": self.segment,
        }


@dataclass
class EventSourcedRepoIsolate:
    """Represents an isolate in a Virtool reference repository."""

    id: UUID
    """The isolate id."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the isolate was not migrated from a legacy repository, this will be `None`.
    """

    name: IsolateName
    """The isolate's source name metadata."""

    sequences: list[EventSourcedRepoSequence]
    """A list of child sequences."""

    @property
    def accessions(self) -> list:
        """Return a list of accessions contained in this isolate"""
        return [sequence.accession for sequence in self.sequences]

    def add_sequence(self, sequence: EventSourcedRepoSequence):
        self.sequences.append(sequence)

    def dict(self):
        return {
            "id": self.id,
            "legacy_id": self.legacy_id,
            "name": self.name.model_dump(),
            "sequences": [sequence.dict() for sequence in self.sequences],
        }


@dataclass
class EventSourcedRepoOTU:
    """Represents an OTU in a Virtool reference repository."""

    id: UUID
    """The OTU id."""

    acronym: str
    """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

    excluded_accessions: list[str]

    isolates: list[EventSourcedRepoIsolate]
    """A list of child isolates."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the OTU was not migrated from a legacy repository, this will be `None`.
    """

    name: str
    """The OTU name (eg. Tobacco mosaic virus)."""

    taxid: int
    """The OTU taxonomy ID."""

    molecule: Molecule | None = None
    """The molecule of this OTU"""

    schema: list | None = None
    """The schema of the OTU"""

    @property
    def accessions(self) -> list:
        """Return a list of accessions contained in this isolate"""
        accessions = []
        for isolate in self.isolates:
            accessions = accessions + isolate.accessions

        return accessions

    @property
    def blocked_accessions(self) -> list:
        return self.accessions + self.excluded_accessions

    def get_isolate(self, isolate_id: UUID) -> EventSourcedRepoIsolate | None:
        """Return the isolate instance associated with a given UUID if it exists,
        else None.
        """
        for isolate in self.isolates:
            if isolate.id == isolate_id:
                return isolate

        return None

    def get_isolate_id(self, type: str, name: str) -> UUID | None:
        """Return an UUID if the name is extant in this OTU."""
        for isolate in self.isolates:
            if isolate.name.type == type and isolate.name.value == name:
                return isolate.id

        return None

    def add_isolate(self, isolate: EventSourcedRepoIsolate):
        self.isolates.append(isolate)

    def dict(self):
        return {
            "id": self.id,
            "acronym": self.acronym,
            "excluded_accessions": self.excluded_accessions,
            "isolates": [isolate.dict() for isolate in self.isolates],
            "legacy_id": self.legacy_id,
            "name": self.name,
            "molecule": self.molecule,
            "schema": self.schema,
            "taxid": self.taxid,
        }
