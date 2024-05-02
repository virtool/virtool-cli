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
from pprint import pprint
from typing import Generator, Type

import arrow
from orjson import orjson
from structlog import get_logger

from virtool_cli.ref.checking import Checker
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
from virtool_cli.ref.utils import DataType, IsolateName, Molecule, pad_zeroes

logger = get_logger("repo")


class EventSourcedRepo:
    def __init__(self, path: Path):
        self.path = path
        """The path to the repo directory."""

        self.last_id = 0

        logger.info("Loading repository")

        for event in self._iter_events():
            if event.id - self.last_id != 1:
                raise ValueError("Event IDs are not sequential.")

            self.last_id = event.id

        self.checker = Checker(self)

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

        src_path = path / "src"
        src_path.mkdir()

        shutil.copytree(
            Path(__file__).parent.parent.parent / "assets/github",
            path / ".github",
        )

        (path / ".cache").mkdir()

        repo_id = uuid.uuid4()

        _write_event(
            src_path,
            1,
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
    def meta(self):
        """The metadata for the repository."""
        for event in self._iter_events():
            if isinstance(event, CreateRepo):
                repo = event.data.model_dump()
                return RepoMeta(**repo, created_at=event.timestamp)

        raise ValueError("No repository creation event found")

    @property
    def src_path(self) -> Path:
        """The path to the repo src directory."""
        return self.path / "src"

    def _iter_events(self):
        for path in sorted(self.src_path.iterdir()):
            yield _read_event_at_path(path)

    def _read_event(self, event_id: int) -> Event:
        return _read_event_at_path(self.src_path / f"{pad_zeroes(event_id)}.json")

    def _write_event(
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

        with open(self.src_path / f"{pad_zeroes(event_id)}.json", "wb") as f:
            f.write(
                orjson.dumps(
                    event.model_dump(by_alias=True),
                    f,
                ),
            )

        self.last_id = event_id

        return event

    def iter_otus(self) -> Generator[EventSourcedRepoOTU, None, None]:
        """Iterate over the OTUs in the repository."""
        otu_event_index = defaultdict(list)

        for event in self._iter_events():
            if type(event) in (
                CreateOTU,
                CreateIsolate,
                CreateSequence,
                ExcludeAccession,
            ):
                otu_event_index[event.query.otu_id].append(event.id)

        for otu_id, _ in otu_event_index.items():
            otu = self.get_otu(otu_id)
            yield otu

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
        self.checker.check_otu_name_exists(name)
        self.checker.check_legacy_id_exists(legacy_id)

        otu_id = uuid.uuid4()

        self._write_event(
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

        return self.get_otu(otu_id)

    def create_isolate(
        self,
        otu_id: uuid.UUID,
        legacy_id: str | None,
        source_name: str,
        source_type: str,
    ):
        isolate_id = uuid.uuid4()

        name = IsolateName(type=source_type, value=source_name)

        self._write_event(
            CreateIsolate,
            CreateIsolateData(id=isolate_id, legacy_id=legacy_id, name=name),
            IsolateQuery(isolate_id=isolate_id, otu_id=otu_id),
        )

        return EventSourcedRepoIsolate(
            id=isolate_id,
            legacy_id=legacy_id,
            sequences=[],
            name=name,
        )

    def create_sequence(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        accession: str,
        definition: str,
        legacy_id: str | None,
        segment: str,
        sequence: str,
    ):
        sequence_id = uuid.uuid4()

        self._write_event(
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

        return EventSourcedRepoSequence(
            id=sequence_id,
            accession=accession,
            definition=definition,
            legacy_id=legacy_id,
            segment=segment,
            sequence=sequence,
        )

    def exclude_accession(self, otu_id: uuid.UUID, accession: str):
        """Exclude an accession for an OTU.

        This accession will not be allowed in the repository in the future.

        :param otu_id: the id of the OTU
        :param accession: the accession to exclude

        """
        self._write_event(
            ExcludeAccession,
            ExcludeAccessionData(accession=accession),
            OTUQuery(otu_id=otu_id),
        )

    def get_otu(self, otu_id: uuid.UUID) -> EventSourcedRepoOTU:
        """Get an OTU by its ID."""
        pprint(list(self._iter_events()))

        event_ids = [
            event.id
            for event in self._iter_events()
            if (
                type(event)
                in (CreateOTU, CreateIsolate, CreateSequence, ExcludeAccession)
                and event.query.otu_id == otu_id
            )
        ]

        first_event_id = event_ids[0]

        event = self._read_event(first_event_id)

        if not isinstance(event, CreateOTU):
            raise ValueError(
                f"The first event ({first_event_id}) for an OTU is not a CreateOTU "
                "event",
            )

        otu = EventSourcedRepoOTU(
            id=event.data.id,
            acronym=event.data.acronym,
            excluded_accessions=[],
            isolates=[],
            legacy_id=event.data.legacy_id,
            molecule=event.data.molecule,
            name=event.data.name,
            schema=event.data.otu_schema,
            taxid=event.data.taxid,
        )

        for event_id in event_ids[1:]:
            event = self._read_event(event_id)

            if isinstance(event, CreateIsolate):
                otu.isolates.append(
                    EventSourcedRepoIsolate(
                        id=event.data.id,
                        legacy_id=event.data.legacy_id,
                        sequences=[],
                        name=event.data.name,
                    ),
                )

            elif isinstance(event, ExcludeAccession):
                otu.excluded_accessions.append(event.data.accession)

            elif isinstance(event, CreateSequence):
                for isolate in otu.isolates:
                    if isolate.id == event.query.isolate_id:
                        isolate.sequences.append(
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
