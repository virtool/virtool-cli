import shutil
from pathlib import Path
from uuid import UUID

from pydantic import BaseModel

from virtool_cli.ref.resources import (
    EventSourcedRepoOTU,
    EventSourcedRepoSequence,
    EventSourcedRepoIsolate,
)
from virtool_cli.ref.snapshot.model import (
    OTUSnapshotOTU,
    OTUSnapshotIsolate,
    OTUSnapshotSequence,
    OTUSnapshotToCIsolate,
    toc_adapter,
)


class OTUSnapshotMeta(BaseModel):
    """Structures metadata about the OTU snapshot itself."""

    at_event: int | None = None
    """The event ID of the last change made to this snapshot."""


class OTUSnapshotToC:
    """Manages the building and loading of the table of contents."""

    def __init__(self, path: Path):
        self.path: Path = path
        """The path to the table of contents."""

    @staticmethod
    def generate_from_otu(
        otu: "EventSourcedRepoOTU",
    ) -> dict[str, OTUSnapshotToCIsolate]:
        """Return a new table of contents from an OTU."""
        toc = {
            f"{isolate.name}": OTUSnapshotToCIsolate(
                id=isolate.id,
                accessions={
                    accession: isolate.get_sequence_by_accession(accession).id
                    for accession in sorted(isolate.accessions)
                },
            )
            for isolate in otu.isolates
        }
        return toc

    def load(self) -> dict[str, OTUSnapshotToCIsolate] | None:
        """Load a table of contents from file."""
        if not self.path.exists():
            return None

        with open(self.path, "rb") as f:
            return toc_adapter.validate_json(f.read())

    def write(self, data: dict[str, OTUSnapshotToCIsolate], indent: int | None = None):
        """Write a table of contents to file."""
        with open(self.path, "wb") as f:
            f.write(toc_adapter.dump_json(data, indent=indent))


class OTUSnapshotDataStore:
    """Stores and retrieves OTU data in snapshot models."""

    def __init__(self, path: Path):
        self.path = path

        if not self.path.exists():
            self.path.mkdir()

    @property
    def contents(self):
        return list(self.path.glob("*.json"))

    def clean(self):
        shutil.rmtree(self.path)
        self.path.mkdir()

    def load_isolate(self, isolate_id: UUID) -> OTUSnapshotIsolate:
        with open(self.path / f"{isolate_id}.json", "rb") as f:
            return OTUSnapshotIsolate.model_validate_json(f.read())

    def cache_isolate(self, isolate: EventSourcedRepoIsolate):
        validated_isolate = OTUSnapshotIsolate(**isolate.dict(exclude_contents=True))
        with open(self.path / f"{isolate.id}.json", "w") as f:
            f.write(validated_isolate.model_dump_json())

    def load_sequence(self, sequence_id: UUID) -> OTUSnapshotSequence:
        with open(self.path / f"{sequence_id}.json", "rb") as f:
            return OTUSnapshotSequence.model_validate_json(f.read())

    def cache_sequence(self, sequence: EventSourcedRepoSequence):
        validated_sequence = OTUSnapshotSequence(**sequence.dict())
        with open(self.path / f"{sequence.id}.json", "w") as f:
            f.write(validated_sequence.model_dump_json())


class OTUSnapshot:
    """Manages snapshot data for a single OTU."""

    def __init__(self, path: Path):
        self.path = path
        """The path of this snapshot's directory."""

        if not self.path.exists():
            self.path.mkdir()

        self._data = OTUSnapshotDataStore(self.path / "data")
        """The data store of this snapshot. Holds isolate and sequence data."""

        self._toc = OTUSnapshotToC(self.path / "toc.json")

        self._metadata_path = self.path / "metadata.json"
        """The path of this snapshot's metadata file."""

        self._metadata = (
            self._load_metadata() if self._metadata_path.exists() else OTUSnapshotMeta()
        )
        """The metadata for this snapshot."""

    @property
    def at_event(self) -> int | None:
        return self._metadata.at_event

    @property
    def _otu_path(self):
        return self.path / "otu.json"

    @property
    def _toc_path(self):
        return self.path / "toc.json"

    def clean(self):
        """Delete and remake OTUSnapshot directory structure."""
        shutil.rmtree(self.path)
        self.path.mkdir()
        self._data.clean()

    def cache(
        self, otu: "EventSourcedRepoOTU", at_event: int | None = None, options=None
    ):
        """Cache an OTU at a given event."""
        self._metadata.at_event = at_event

        self._cache_otu(otu)

        for isolate in otu.isolates:
            self._data.cache_isolate(isolate)

            for sequence in isolate.sequences:
                self._data.cache_sequence(sequence)

        self._toc.write(data=OTUSnapshotToC.generate_from_otu(otu))

        self._write_metadata(indent=None)

    def _cache_otu(self, otu: EventSourcedRepoOTU):
        validated_otu = OTUSnapshotOTU(**otu.dict(exclude_contents=True))
        with open(self._otu_path, "w") as f:
            f.write(validated_otu.model_dump_json())

    def cache_isolate(self, isolate, at_event: int | None = None, options=None):
        self._metadata.at_event = at_event

        toc = self._toc.load()
        toc[f"{isolate.name}"] = OTUSnapshotToCIsolate(
            id=isolate.id,
            accessions={
                accession: isolate.get_sequence_by_accession(accession).id
                for accession in sorted(isolate.accessions)
            },
        )
        self._toc.write(data=toc)

        self._data.cache_isolate(isolate)

        self._write_metadata(indent=None)

    def cache_sequence(
        self, sequence, isolate_id: UUID, at_event: int | None = None, options=None
    ):
        """Add a new sequence to the Snapshot."""
        self._metadata.at_event = at_event

        toc = self._toc.load()
        for key in toc:
            if toc[key].id == isolate_id:
                toc[key].accessions[sequence.accession] = sequence.id
                break
        self._toc.write(toc)

        self._data.cache_sequence(sequence)

        self._write_metadata(indent=None)

    def load(self) -> "EventSourcedRepoOTU":
        """Load an OTU from the snapshot."""
        with open(self._otu_path, "rb") as f:
            otu_structure = OTUSnapshotOTU.model_validate_json(f.read())

        toc = self._toc.load()

        isolates = []
        for key in toc:
            isolate_entry = toc[key]

            isolate_structure = self._data.load_isolate(isolate_entry.id)

            sequences = []

            for accession in toc[key].accessions:
                sequence_id = toc[key].accessions[accession]
                sequence_structure = self._data.load_sequence(sequence_id)

                sequence = EventSourcedRepoSequence(**sequence_structure.model_dump())

                sequences.append(sequence)

            isolate_dict = isolate_structure.model_dump()
            isolate_dict["uuid"] = isolate_dict.pop("id")

            isolate = EventSourcedRepoIsolate(**isolate_dict, sequences=sequences)

            isolates.append(isolate)

        otu_dict = otu_structure.model_dump(by_alias=True)
        otu_dict["uuid"] = otu_dict.pop("id")

        return EventSourcedRepoOTU(**otu_dict, isolates=isolates)

    def _write_metadata(self, indent) -> None:
        with open(self._metadata_path, "w") as f:
            f.write(self._metadata.model_dump_json(indent=indent))

    def _load_metadata(self) -> OTUSnapshotMeta | None:
        if self._metadata_path.exists():
            with open(self._metadata_path, "rb") as f:
                metadata = OTUSnapshotMeta.model_validate_json(f.read())

            return metadata

        return None
