import datetime
import dataclasses
from dataclasses import dataclass
from pydantic import BaseModel
from uuid import UUID

from virtool_cli.ref.utils import DataType, IsolateName, IsolateNameKey, Molecule


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


class EventSourcedRepoIsolate:
    """Represents an isolate in a Virtool reference repository."""

    def __init__(
        self,
        id: UUID,
        name: IsolateName,
        legacy_id: str | None = None,
        _sequences_by_accession: dict | None = None,
    ):
        self.id = id
        """The isolate id."""

        self.name = name
        """The isolate's source name metadata."""

        self._sequences_by_accession = (
            {} if _sequences_by_accession is None else _sequences_by_accession
        )
        """A dictionary of sequences indexed by accession"""

        self.legacy_id = legacy_id
        """A string based ID carried over from a legacy Virtool reference repository.
    
        It the isolate was not migrated from a legacy repository, this will be `None`.
        """

    @property
    def sequences(self) -> list[EventSourcedRepoSequence]:
        """Return a list of child sequences."""
        return list(self._sequences_by_accession.values())

    @property
    def accessions(self) -> set[str]:
        """Return a set of accessions contained in this isolate."""
        return set(self._sequences_by_accession.keys())

    def __repr__(self) -> str:
        """Return a shorthand representation of the isolate's contents"""
        return (
            "EventSourcedRepoIsolate("
            + f"{self.id}, type={self.name.type}, name={self.name.value}, "
            + f"accessions={self.accessions})"
        )

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


class EventSourcedRepoOTU:
    """Represents an OTU in a Virtool reference repository."""

    def __init__(
        self,
        id: UUID,
        taxid: int,
        name: str,
        acronym: str = "",
        molecule: Molecule | None = None,
        legacy_id: str | None = None,
        schema: list | None = None,
        excluded_accessions: list | None = None,
        _isolates_by_id: dict | None = None,
    ):
        self.id = id
        """The OTU id."""

        self.taxid = taxid
        """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

        self.name = name
        """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

        self.acronym = acronym
        """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

        self.legacy_id = legacy_id
        """A string based ID carried over from a legacy Virtool reference repository."""

        self.molecule = molecule
        """The molecule of this OTU"""

        self.schema = schema
        """The schema of the OTU"""

        self.excluded_accessions = (
            set() if excluded_accessions is None else set(excluded_accessions)
        )
        """A set of accessions that should not be retrieved in future fetch operations"""

        self._isolates_by_id = {} if _isolates_by_id is None else _isolates_by_id
        """A dictionary of isolates indexed by isolate UUID"""

    @property
    def isolates(self) -> list:
        """Returns a list of child isolates"""
        return list(self._isolates_by_id.values())

    @property
    def accessions(self) -> set:
        """Return a set of accessions contained in this isolate"""
        accessions = set()
        for isolate in self.isolates:
            accessions.update(isolate.accessions)

        return accessions

    def __repr__(self) -> str:
        """Return a shorthand representation of the OTU's contents"""
        return (
            "EventSourcedRepoOTU("
            + f"{self.id}, taxid={self.taxid}, name={self.name}, "
            + f"accessions={self.accessions})"
        )

    @property
    def blocked_accessions(self) -> set:
        """Returns a set of accessions to be blocked from fetches,
        i.e. accessions that have already been added and excluded accessions."""
        return self.accessions.union(self.excluded_accessions)

    def add_isolate(self, isolate: EventSourcedRepoIsolate):
        self._isolates_by_id[isolate.id] = isolate

    def get_isolate(self, isolate_id: UUID) -> EventSourcedRepoIsolate | None:
        """Return the isolate instance associated with a given UUID if it exists,
        else None.
        """
        if isolate_id in self._isolates_by_id:
            return self._isolates_by_id[isolate_id]

        return None

    def get_isolate_id_by_name(self, name: IsolateName | IsolateNameKey) -> UUID | None:
        """Return an UUID if the name is extant in this OTU."""
        for isolate in self.isolates:
            if isolate.name == name:
                return isolate.id

        return None

    def dict(self):
        """Return data in JSON-ready form"""
        return {
            "id": self.id,
            "acronym": self.acronym,
            "excluded_accessions": list(self.excluded_accessions),
            "isolates": [isolate.dict() for isolate in self.isolates],
            "legacy_id": self.legacy_id,
            "name": self.name,
            "molecule": self.molecule,
            "schema": self.schema,
            "taxid": self.taxid,
        }
