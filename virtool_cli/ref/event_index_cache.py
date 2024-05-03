import shutil
import dataclasses
from pathlib import Path
from uuid import UUID
from typing import Annotated
from collections.abc import Generator

import orjson
from pydantic import Field, ValidationError
from pydantic.dataclasses import dataclass

from structlog import get_logger

logger = get_logger("repo.cache.event_index")


@dataclass
class OTUEventCache:
    otu_id: UUID

    at_event: Annotated[int, Field(ge=1)]
    """The latest event ID when this cache index was last updated"""

    events: list[int]
    """A list of event IDs"""


class EventIndexCacheError(Exception):
    """Event index cache is corrupted."""


class EventIndexCache:
    """Maintains an index of event records for each OTU UUID in the repo"""

    def __init__(self, path: Path):
        self.path = path

        self.path.mkdir(exist_ok=True)

    def iter_event_records(self) -> Generator[OTUEventCache, None, None]:
        """Iterates through each OTU in the index"""
        for subpath in self.path.iterdir():
            if subpath.suffix == ".json":
                otu_id = UUID(subpath.stem)
                yield self.load_otu_events(otu_id)

    def list_otu_ids(self) -> list[UUID]:
        """Returns a list of OTU Ids with extant data in the cache."""
        return [
            UUID(subpath.stem)
            for subpath in self.path.iterdir()
            if subpath.suffix == ".json"
        ]

    def cache_index(self, event_index: dict[UUID, list[int]], last_id: int):
        for otu_id in event_index:
            self.cache_otu_events(otu_id, event_index[otu_id], last_id=last_id)

    def clear(self):
        """Clear and reset cache."""
        shutil.rmtree(self.path)
        self.path.mkdir()

    def cache_otu_events(self, otu_id: UUID, event_list: list[int], last_id: int = 1):
        """
        :param otu_id: A Virtool OTU Id
        :param event_list: A list of event IDs
        :param last_id: The last added Id in the event store at point of caching
        """
        index_data = OTUEventCache(otu_id=otu_id, at_event=last_id, events=event_list)

        try:
            with open(self.path / f"{otu_id}.json", "wb") as f:
                f.write(orjson.dumps(dataclasses.asdict(index_data)))
        finally:
            logger.debug(
                "OTU event index updated", otu_id=str(otu_id), events=event_list
            )

    def load_otu_events(self, otu_id: UUID) -> OTUEventCache | None:
        """
        Takes a requested OTU Id and returns an OTUEventCache(otu_id, at_event, events)
        if possible, else None if the event record does not exist.

        If the cached data is not compatible with the repo's event store,
        raise an EventIndexCacheError.

        :param otu_id: A Virtool OTU Id
        :return: OTUEventIndex(at_event, events) if a valid event list
            can be found in the cache, else None
        """
        cached_path = self.path / f"{otu_id}.json"

        if not cached_path.exists():
            return None

        with open(cached_path, "rb") as f:
            cached_data = orjson.loads(f.read())

        try:
            otu_index_cache = OTUEventCache(**cached_data)
        except ValidationError:
            raise EventIndexCacheError(
                "Bad Index: Events could not be retrieved from cache"
            )

        return otu_index_cache

    def clear_cached_otu_events(self, otu_id: UUID):
        """Delete a given OTU's cached data from the cache."""
        (self.path / f"{otu_id}.json").unlink(missing_ok=True)
