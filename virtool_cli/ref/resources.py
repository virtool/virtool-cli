import datetime
from dataclasses import dataclass
from uuid import UUID

from virtool_cli.ref.utils import DataType, IsolateName, Molecule


@dataclass
class RepoMeta:
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

    molecule: Molecule
    """The molecule of this OTU"""

    taxid: int
    """The OTU taxonomy ID."""

    schema: list

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
