from dataclasses import dataclass
from pathlib import Path
from uuid import UUID

import orjson
from structlog import get_logger

logger = get_logger("repo.cache.event_index")


@dataclass
class EventIndexItem:
    at_event: int
    """The latest event ID when this cache index was last updated"""

    event_ids: list[int]
    """A list of event IDs"""

    otu_id: UUID
    """The OTU ID this cache index is associated with."""


class EventIndexError(Exception):
    """Event index cache is corrupted."""


class EventIndex:
    """Maintains an index of event records for each OTU UUID in the repo"""

    def __init__(self, path: Path):
        self.path = path
        self.path.mkdir(exist_ok=True)

    def set(self, otu_id: UUID, event_ids: list[int], last_id: int):
        """Cache the event list for a given OTU.

        :param otu_id: an OTU ID
        :param event_ids: the list of event IDs
        :param last_id: the ID of the event last added to the event store
        """
        if last_id < 1:
            raise ValueError("last_id must be greater than 0")

        with open(self.path / f"{otu_id}.json", "wb") as f:
            f.write(
                orjson.dumps(
                    {
                        "at_event": last_id,
                        "event_ids": sorted(event_ids),
                    },
                ),
            )

    def load(self):
        """Cache an event index dictionary"""
        event_index = {}

        for path in self.path.iterdir():
            if path.suffix == ".json":
                otu_id = UUID(path.stem)
                event_index[otu_id] = self.get(otu_id).event_ids

        return event_index

    def get(self, otu_id: UUID) -> EventIndexItem | None:
        """Takes a requested OTU Id and returns an OTUEventCache(otu_id, at_event, events)
        if possible, else None if the event record does not exist.

        If the cached data is not compatible with the repo's event store,
        raise an EventIndexCacheError.

        :param otu_id: A Virtool OTU Id
        :return: OTUEventIndex(at_event, events) if a valid event list
            can be found in the cache, else None
        """
        path = self.path / f"{otu_id}.json"

        try:
            with open(path, "rb") as f:
                return EventIndexItem(
                    **{**orjson.loads(f.read()), "otu_id": otu_id},
                )
        except FileNotFoundError:
            return None
