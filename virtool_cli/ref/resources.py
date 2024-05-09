import datetime
import dataclasses
from dataclasses import dataclass
from pydantic import BaseModel
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

    _sequences_by_accession: dict = dataclasses.field(default_factory=dict)
    """A dictionary of sequences indexed by accession"""

    @property
    def sequences(self) -> list[EventSourcedRepoSequence]:
        """Return a list of child sequences."""
        return list(self._sequences_by_accession.values())

    @property
    def accession_set(self) -> set[str]:
        """Return a set of accessions contained in this isolate."""
        return set(self._sequences_by_accession.keys())

    def add_sequence(self, sequence: EventSourcedRepoSequence):
        """Add a new sequence to private dictionary"""
        self._sequences_by_accession[sequence.accession] = sequence

    def get_sequence_by_accession(
        self, accession: str
    ) -> EventSourcedRepoSequence | None:
        """Return a sequence with the given accession if it exists in the isolate,
        else None"""
        if accession in self._sequences_by_accession:
            return self._sequences_by_accession[accession]

        return None

    def dict(self) -> dict:
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

    excluded_accessions: set[str]

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

    _isolates_by_id: dict = dataclasses.field(default_factory=dict)
    """A dictionary of isolates indexed by isolate UUID"""

    @property
    def isolates(self) -> list:
        """Returns a list of child isolates"""
        return list(self._isolates_by_id.values())

    @property
    def accession_set(self) -> set:
        """Return a set of accessions contained in this isolate"""
        accessions = set()
        for isolate in self.isolates:
            accessions.update(isolate.accession_set)

        return accessions

    @property
    def blocked_accession_set(self) -> set:
        """Returns a set of accessions to be blocked from fetches,
        i.e. accessions that have already been added and excluded accessions."""
        return self.accession_set.union(self.excluded_accessions)

    def add_isolate(self, isolate: EventSourcedRepoIsolate):
        if isolate.id in self._isolates_by_id:
            raise ValueError("Isolate already exists")

        self._isolates_by_id[isolate.id] = isolate

    def get_isolate(self, isolate_id: UUID) -> EventSourcedRepoIsolate | None:
        """Return the isolate instance associated with a given UUID if it exists,
        else None.
        """
        if isolate_id in self._isolates_by_id:
            return self._isolates_by_id[isolate_id]

        return None

    def get_isolate_id_by_name(self, name: IsolateName) -> UUID | None:
        """Return an UUID if the name is extant in this OTU."""
        for isolate in self.isolates:
            if isolate.name == name:
                return isolate.id

        return None

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
