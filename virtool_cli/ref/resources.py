import datetime
from dataclasses import dataclass
from uuid import UUID

from virtool_cli.ref.utils import DataType


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

    sequence: str
    """The sequence."""

    segment: str
    """The sequence segment."""

    def dict(self):
        return {
            "id": self.id,
            "accession": self.accession,
            "definition": self.definition,
            "sequence": self.sequence,
            "segment": self.segment,
        }


@dataclass
class EventSourcedRepoIsolate:
    """Represents an isolate in a Virtool reference repository."""

    id: UUID
    """The isolate id."""

    source_name: str
    """The isolate's source name."""

    source_type: str
    """The isolate's source type."""

    sequences: list[EventSourcedRepoSequence]
    """A list of child sequences."""

    def add_sequence(self, sequence: EventSourcedRepoSequence):
        self.sequences.append(sequence)

    def dict(self):
        return {
            "id": self.id,
            "source_name": self.source_name,
            "source_type": self.source_type,
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

    name: str
    """The OTU name (eg. Tobacco mosaic virus)."""

    schema: list

    taxid: int
    """The OTU taxonomy ID."""

    def add_isolate(self, isolate: EventSourcedRepoIsolate):
        self.isolates.append(isolate)

    def dict(self):
        return {
            "id": self.id,
            "acronym": self.acronym,
            "excluded_accessions": self.excluded_accessions,
            "isolates": [isolate.dict() for isolate in self.isolates],
            "name": self.name,
            "schema": self.schema,
            "taxid": self.taxid,
        }
