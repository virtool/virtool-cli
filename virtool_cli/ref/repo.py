"""An access layer for a Virtool event-sourced reference repository.

This is a work in progress.

**Check that OTUs have only one representative (default) isolate.**

The default isolate is set using `rep_isolate` in the `CreateOTU` event. Therefore only
one isolate in an OTU can be default.

TODO: Check if excluded accessions exist in the repo.
TODO: Check for accessions filed in wrong isolates.
TODO: Check for accession conflicts.

"""

import shutil
import uuid
from collections import defaultdict
from pathlib import Path
from typing import Generator, Type

import arrow
from orjson import orjson
from structlog import get_logger

from virtool_cli.ref.events import (
    CreateIsolate,
    CreateIsolateData,
    CreateOTU,
    CreateOTUData,
    CreateRepo,
    CreateRepoData,
    CreateSequence,
    CreateSequenceData,
    Event,
    EventData,
    EventQuery,
    ExcludeAccession,
    ExcludeAccessionData,
    IsolateQuery,
    OTUQuery,
    RepoQuery,
    SequenceQuery,
)
from virtool_cli.ref.resources import (
    EventSourcedRepoIsolate,
    EventSourcedRepoOTU,
    EventSourcedRepoSequence,
    RepoMeta,
)
from virtool_cli.ref.snapshot.index import SnapshotIndex
from virtool_cli.ref.event_index_cache import EventIndexCache, EventIndexCacheError
from virtool_cli.utils.models import Molecule
from virtool_cli.ref.utils import DataType, IsolateName, pad_zeroes

logger = get_logger("repo")


OTU_EVENT_TYPES = (CreateOTU, CreateIsolate, CreateSequence, ExcludeAccession)


class EventSourcedRepo:
    def __init__(self, path: Path):
        self.path = path
        """The path to the repo directory."""

        logger.info("Loading repository")

        self.cache_path = self.path / ".cache"
        """The path to the cache subdirectory."""

        self._event_store = EventStore(self.path / "src")
        """The event store of the event sourced repository."""

        self._snapshot_path = path / ".cache/snapshot"
        """The path to the snapshot cache subdirectory."""

        self._event_index_cache = EventIndexCache(self.cache_path / "event_index")
        """The event index cache of the event sourced repository."""

        self._snapshotter = (
            SnapshotIndex(path=self._snapshot_path)
            if self._snapshot_path.exists()
            else SnapshotIndex.new(path=self._snapshot_path, metadata=self.meta)
        )
        """The snapshot index. Maintains and caches the read model of the Repo."""

        # Take a new snapshot if no existing data is found.
        if not self._snapshotter.otu_ids:
            logger.debug("No snapshot data found. Building new snapshot...")
            self.snapshot()

        logger.info("Finished loading repository", event_count=self.last_id)

    @classmethod
    def new(cls, data_type: DataType, name: str, path: Path, organism: str):
        """Create a new reference repository."""
        if path.is_file():
            raise ValueError("The target path is a file")

        path.mkdir(parents=True, exist_ok=True)

        if any(path.iterdir()):
            raise ValueError("The target path is not empty")

        with open(path / ".gitignore", "w") as f:
            f.write(".cache\n")

        shutil.copytree(
            Path(__file__).parent.parent.parent / "assets/github",
            path / ".github",
        )

        (path / ".cache").mkdir()

        repo_id = uuid.uuid4()

        _src = EventStore(path / "src")
        _src.write_event(
            CreateRepo,
            CreateRepoData(
                id=repo_id,
                data_type=data_type,
                name=name,
                organism=organism,
            ),
            RepoQuery(repository_id=repo_id),
        )

        return EventSourcedRepo(path)

    @property
    def last_id(self):
        """The id of the most recently added event in the event store."""
        return self._event_store.last_id

    @property
    def meta(self):
        """The metadata for the repository."""
        for event in self._event_store.iter_events():
            if isinstance(event, CreateRepo):
                repo = event.data.model_dump()
                return RepoMeta(**repo, created_at=event.timestamp)

        raise ValueError("No repository creation event found")

    @property
    def src_path(self) -> Path:
        """The path to the repo src directory."""
        return self._event_store.path

    @property
    def taxids(self) -> set:
        """Extant Taxonomy ids in the read model"""
        return self._snapshotter.taxids

    @property
    def accessions(self) -> set:
        """Extant accessions in the read model"""
        return self._snapshotter.accessions

    def _get_event_index(self) -> dict[uuid.UUID, list[int]]:
        """Get the current event index from the event store,
        binned and indexed by OTU Id."""
        otu_event_index = defaultdict(list)

        for event in self._event_store.iter_events():
            if type(event) in OTU_EVENT_TYPES:
                otu_event_index[event.query.otu_id].append(event.id)

        return otu_event_index

    def _get_event_index_after_start(
        self, start: int = 1
    ) -> dict[uuid.UUID, list[int]]:
        """Get the current event index, binned and indexed by OTU ID"""
        otu_event_index = defaultdict(list)

        for event in self._event_store.iter_events_from_index(start):
            if type(event) in OTU_EVENT_TYPES:
                otu_event_index[event.query.otu_id].append(event.id)

        return otu_event_index

    def snapshot(self):
        """Create a snapshot using all the OTUs in the event store."""
        self._snapshotter.snapshot(
            self.get_all_otus(ignore_cache=True),
            at_event=self.last_id,
            indent=True,
        )

    def iter_otus(self) -> Generator[EventSourcedRepoOTU, None, None]:
        """Iterate over the OTUs in the snapshot."""
        for otu_id in self._snapshotter.otu_ids:
            otu = self._snapshotter.get_otu(otu_id)
            yield otu

    def get_all_otus(self, ignore_cache: bool = False) -> list[EventSourcedRepoOTU]:
        """Retrieve all OTUs from the event store and return as a list."""
        if ignore_cache:
            event_index = self._get_event_index()
        else:
            event_index = self._event_index_cache.load_index()

        otus = []
        for otu_id in event_index:
            if otu := self.get_otu(otu_id, ignore_cache):
                otus.append(otu)

        return otus

    def create_otu(
        self,
        acronym: str,
        legacy_id: str | None,
        name: str,
        molecule: Molecule | None,
        schema: [],
        taxid: int,
    ):
        """Create an OTU."""
        if taxid in self._snapshotter.taxids:
            raise ValueError(
                f"OTU already exists as {self._snapshotter.index_by_taxid[taxid]}",
            )
        if name in self._snapshotter.index_by_name:
            raise ValueError(f"An OTU with the name '{name}' already exists")
        if legacy_id in self._snapshotter.index_by_legacy_id:
            raise ValueError(f"An OTU with the legacy ID '{legacy_id}' already exists")

        otu_logger = logger.bind(taxid=taxid, name=name, legacy_id=legacy_id)
        otu_logger.info(f"Creating new OTU for Taxonomy ID {taxid}...")

        otu_id = uuid.uuid4()

        event = self._event_store.write_event(
            CreateOTU,
            CreateOTUData(
                id=otu_id,
                acronym=acronym,
                legacy_id=legacy_id,
                name=name,
                molecule=molecule,
                schema=schema,
                rep_isolate=None,
                taxid=taxid,
            ),
            OTUQuery(otu_id=otu_id),
        )

        otu_logger.debug("OTU written", event_id=event.id, otu_id=str(otu_id))

        otu = self.get_otu(otu_id)

        self._snapshotter.cache_otu(otu, at_event=self.last_id)

        return otu

    def create_isolate(
        self,
        otu_id: uuid.UUID,
        legacy_id: str | None,
        source_name: str,
        source_type: str,
    ) -> EventSourcedRepoIsolate | None:
        """Create and return a new isolate within the given OTU.
        If the isolate name already exists, return None."""
        otu = self.get_otu(otu_id, ignore_cache=False)

        name = IsolateName(**{"type": source_type, "value": source_name})
        if otu.get_isolate_id_by_name(name) is not None:
            logger.warning(
                "An isolate by this name already exists", isolate_name=str(name)
            )
            return None

        isolate_id = uuid.uuid4()

        event = self._event_store.write_event(
            CreateIsolate,
            CreateIsolateData(id=isolate_id, legacy_id=legacy_id, name=name),
            IsolateQuery(isolate_id=isolate_id, otu_id=otu_id),
        )

        logger.debug(
            "Isolate written",
            event_id=event.id,
            isolate_id=str(isolate_id),
            name=str(name),
        )

        isolate = EventSourcedRepoIsolate(
            uuid=isolate_id,
            legacy_id=legacy_id,
            name=name,
            sequences=[],
        )

        otu.add_isolate(isolate)

        self._snapshotter.cache_otu(otu, at_event=self.last_id)

        return isolate

    def create_sequence(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        accession: str,
        definition: str,
        legacy_id: str | None,
        segment: str,
        sequence: str,
    ) -> EventSourcedRepoSequence | None:
        """Create and return a new sequence within the given OTU.
        If the accession already exists in this OTU, return None."""
        otu = self.get_otu(otu_id, ignore_cache=False)

        if accession in otu.accessions:
            logger.warning(
                "This accession already exists in the OTU.",
                accession=accession,
                otu_id=str(otu_id),
            )
            return None

        sequence_id = uuid.uuid4()

        event = self._event_store.write_event(
            CreateSequence,
            CreateSequenceData(
                id=sequence_id,
                accession=accession,
                definition=definition,
                legacy_id=legacy_id,
                segment=segment,
                sequence=sequence,
            ),
            SequenceQuery(
                otu_id=otu_id,
                isolate_id=isolate_id,
                sequence_id=sequence_id,
            ),
        )

        logger.debug(
            "Sequence written",
            event_id=event.id,
            sequence_id=str(sequence_id),
            accession=accession,
        )

        sequence = EventSourcedRepoSequence(
            id=sequence_id,
            accession=accession,
            definition=definition,
            legacy_id=legacy_id,
            segment=segment,
            sequence=sequence,
        )

        otu.add_sequence(sequence, isolate_id)

        self._snapshotter.cache_otu(otu, at_event=self.last_id)

        return sequence

    def exclude_accession(self, otu_id: uuid.UUID, accession: str):
        """Exclude an accession for an OTU.

        This accession will not be allowed in the repository in the future.

        :param otu_id: the id of the OTU
        :param accession: the accession to exclude

        """
        self._event_store.write_event(
            ExcludeAccession,
            ExcludeAccessionData(accession=accession),
            OTUQuery(otu_id=otu_id),
        )

    def read_otu(self, otu_id: uuid.UUID) -> EventSourcedRepoOTU | None:
        """Return an OTU corresponding to a UUID if found in snapshot, else None"""
        logger.debug("Loading OTU from snapshot...", otu_id=str(otu_id))

        return self._snapshotter.load_otu(otu_id)

    def read_otu_by_taxid(self, taxid: int) -> EventSourcedRepoOTU | None:
        """Return an OTU corresponding to a Taxonomy ID if found in snapshot, else None"""
        logger.debug("Loading OTU from snapshot...", taxid=taxid)

        return self._snapshotter.load_otu_by_taxid(taxid)

    def get_otu(
        self, otu_id: uuid.UUID, ignore_cache: bool = False
    ) -> EventSourcedRepoOTU | None:
        """Return an OTU corresponding with a given OTU Id if it exists, else None."""
        logger.debug("Getting OTU from events...", otu_id=str(otu_id))
        event_ids = self._get_otu_events(otu_id, ignore_cache)

        if event_ids:
            return self._rehydrate_otu(event_ids)

        return None

    def get_otu_by_taxid(
        self, taxid: int, ignore_cache: bool = False
    ) -> EventSourcedRepoOTU | None:
        """Return an OTU corresponding with a given OTU Id if it exists, else None"""
        if (otu_id := self._snapshotter.index_by_taxid.get(taxid)) is not None:
            otu = self.get_otu(otu_id, ignore_cache=ignore_cache)
            return otu

        return None

    def _rehydrate_otu(self, event_ids: list[int]) -> EventSourcedRepoOTU:
        """Rebuilds OTU data from a list of event IDs"""
        first_event_id = sorted(event_ids)[0]

        event = self._event_store.read_event(first_event_id)

        if not isinstance(event, CreateOTU):
            raise ValueError(
                f"The first event ({first_event_id}) for an OTU is not a CreateOTU "
                "event",
            )

        otu = EventSourcedRepoOTU(
            uuid=event.data.id,
            acronym=event.data.acronym,
            excluded_accessions=[],
            legacy_id=event.data.legacy_id,
            name=event.data.name,
            taxid=event.data.taxid,
            molecule=event.data.molecule,
            schema=event.data.otu_schema,
        )

        for event_id in event_ids[1:]:
            event = self._event_store.read_event(event_id)

            if isinstance(event, CreateIsolate):
                otu.add_isolate(
                    EventSourcedRepoIsolate(
                        uuid=event.data.id,
                        legacy_id=event.data.legacy_id,
                        name=event.data.name,
                    ),
                )

            elif isinstance(event, ExcludeAccession):
                otu.excluded_accessions.add(event.data.accession)

            elif isinstance(event, CreateSequence):
                for isolate in otu.isolates:
                    if isolate.id == event.query.isolate_id:
                        isolate.add_sequence(
                            EventSourcedRepoSequence(
                                id=event.data.id,
                                accession=event.data.accession,
                                definition=event.data.definition,
                                legacy_id=event.data.legacy_id,
                                segment=event.data.segment,
                                sequence=event.data.sequence,
                            ),
                        )

        return otu

    def _get_otu_metadata(self, event_ids: list[int]) -> dict | None:
        """Retrieves OTU metadata from a list of event IDs"""
        if not event_ids:
            return None
        event_ids.sort()
        first_event_id = event_ids[0]

        event = self._event_store.read_event(first_event_id)

        if not isinstance(event, CreateOTU):
            raise ValueError(
                f"The first event ({first_event_id}) for an OTU is not a CreateOTU "
                "event",
            )

        return {
            "id": event.data.id,
            "acronym": event.data.acronym,
            "legacy_id": event.data.legacy_id,
            "name": event.data.name,
            "taxid": event.data.taxid,
        }

    def _get_otu_events(
        self, otu_id: uuid.UUID, ignore_cache: bool = False
    ) -> list[int]:
        """
        Returns an up-to-date list of events associated with this OTU Id.

        If ignore_cache, loads the OTU's event index cache and makes sure
        the results are up to date before returning the list.

        If cache load fails or ignore_cache is False, generates a new event list.
        """
        otu_logger = logger.bind(otu_id=str(otu_id), ignore_cache=ignore_cache)

        if not ignore_cache:
            try:
                otu_event_list = self._load_otu_events_from_cache_and_update(otu_id)

                if otu_event_list:
                    otu_logger.debug(
                        "Cached events found",
                        event_ids=otu_event_list,
                        last_event=self.last_id,
                    )
                    return otu_event_list

                otu_logger.debug(
                    "Event index cache was empty.",
                    event_ids=otu_event_list,
                    last_event=self.last_id,
                )

            except EventIndexCacheError as e:
                logger.warning(e, last_event=self.last_id)

                logger.debug("Deleting bad index...")
                self._event_index_cache.clear_cached_otu_events(otu_id)

        otu_logger.debug("Searching event store for matching events...")

        event_ids = [
            event.id
            for event in self._event_store.iter_events()
            if (type(event) in OTU_EVENT_TYPES and event.query.otu_id == otu_id)
        ]

        otu_logger.debug("Writing events to cache...", events=event_ids)

        self._event_index_cache.cache_otu_events(
            otu_id, event_ids, last_id=self.last_id
        )

        return event_ids

    def _load_otu_events_from_cache_and_update(self, otu_id: uuid.UUID) -> list[int]:
        """Gets OTU events from the event index cache,
        updates the list it is not up to date and returns the list.
        """
        otu_logger = logger.bind(otu_id=str(otu_id))

        cached_otu_index = None
        try:
            cached_otu_index = self._event_index_cache.load_otu_events(otu_id)

        except EventIndexCacheError as e:
            logger.error(e)

            logger.warning("Deleting bad index...")

            self._event_index_cache.clear_cached_otu_events(otu_id)

        if cached_otu_index is not None:
            if cached_otu_index.at_event == self.last_id:
                return cached_otu_index.events

            if cached_otu_index.at_event > self.last_id:
                raise EventIndexCacheError(
                    "Bad Index: "
                    + "Cached event index is greater than current repo's last ID"
                )

            # Update event list
            otu_logger.debug(
                "Cached event list is out of date",
                cached_events=cached_otu_index.events,
                cache_at=cached_otu_index.at_event,
                last_event=self.last_id,
            )

            otu_event_list = cached_otu_index.events

            for event in self._event_store.iter_events_from_index(
                start=cached_otu_index.at_event
            ):
                if type(event) in OTU_EVENT_TYPES and event.id not in otu_event_list:
                    otu_event_list.append(event.id)

            otu_event_list.sort()

            otu_logger.debug(
                "Added new events to event list.",
                cached_events=cached_otu_index.events,
                updated_events=otu_event_list,
            )

            self._event_index_cache.cache_otu_events(
                otu_id, otu_event_list, last_id=self.last_id
            )

            return otu_event_list

        return []


class EventStore:
    """Interface for the event store"""

    def __init__(self, path: Path):
        path.mkdir(exist_ok=True)

        self.path = path
        """The path to the event store directory."""

        self.last_id = 0
        """The id of the latest event."""

        # Check that all events are present and set .last_id to the latest event.
        for event_id in self.event_ids:
            if event_id - self.last_id != 1:
                raise ValueError("Event IDs are not sequential.")

            self.last_id = event_id

    @property
    def event_ids(self) -> list:
        event_ids = []
        for event_path in self.path.iterdir():
            try:
                event_ids.append(int(event_path.stem))
            except ValueError:
                continue
        event_ids.sort()

        return event_ids

    def iter_events(self, reverse: bool = False):
        for path in sorted(self.path.glob("*.json"), reverse=reverse):
            if path.stem == "meta":
                continue
            yield EventStore._read_event_at_path(path)

    def iter_events_from_index(self, start: int = 1) -> Generator[Event, None, None]:
        """Iterates through events in the src directory using the index id.

        :param start: The event ID to be read first.
            Setting start > 1 will begin the iterator
            from the middle of the event store.
        """
        if start < 1:
            raise IndexError("Start index cannot be <1")

        if start > self.last_id:
            raise IndexError(f"Start index cannot be >{self.last_id}")

        for event_index in range(start, self.last_id + 1):
            yield self.read_event(event_index)

    def read_event(self, event_id: int) -> Event:
        return EventStore._read_event_at_path(
            self.path / f"{pad_zeroes(event_id)}.json"
        )

    def write_event(
        self,
        cls: Type[Event],
        data: EventData,
        query: EventQuery,
    ) -> Event:
        """Write a new event to the repository."""
        event_id = self.last_id + 1

        event = cls(
            id=event_id,
            data=data,
            query=query,
            timestamp=arrow.utcnow().naive,
        )

        with open(self.path / f"{pad_zeroes(event_id)}.json", "wb") as f:
            f.write(
                orjson.dumps(
                    event.model_dump(by_alias=True),
                    f,
                ),
            )

        self.last_id = event_id

        return event

    @staticmethod
    def _read_event_at_path(path: Path) -> Event:
        with open(path, "rb") as f:
            loaded = orjson.loads(f.read())

            match loaded["type"]:
                case "CreateRepo":
                    return CreateRepo(**loaded)
                case "CreateOTU":
                    return CreateOTU(**loaded)
                case "CreateIsolate":
                    return CreateIsolate(**loaded)
                case "CreateSequence":
                    return CreateSequence(**loaded)
                case "ExcludeAccession":
                    return ExcludeAccession(**loaded)

            raise ValueError(f"Unknown event type: {loaded['type']}")


def _write_event(
    src_path: Path,
    event_id: int,
    cls: Type[Event],
    data: EventData,
    query: EventQuery,
):
    event = cls(
        id=event_id,
        data=data,
        query=query,
        timestamp=arrow.utcnow().naive,
    )

    with open(src_path / f"{pad_zeroes(event_id)}.json", "wb") as f:
        f.write(
            orjson.dumps(
                event.model_dump(),
                f,
            ),
        )

    return event


def _read_event_at_path(path: Path) -> Event:
    with open(path, "rb") as f:
        loaded = orjson.loads(f.read())

        match loaded["type"]:
            case "CreateRepo":
                return CreateRepo(**loaded)
            case "CreateOTU":
                return CreateOTU(**loaded)
            case "CreateIsolate":
                return CreateIsolate(**loaded)
            case "CreateSequence":
                return CreateSequence(**loaded)
            case "ExcludeAccession":
                return ExcludeAccession(**loaded)

        raise ValueError(f"Unknown event type: {loaded['type']}")
