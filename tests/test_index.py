"""Test for ``EventIndex`` class."""

import uuid
from pathlib import Path

from virtool_cli.ref.event_index_cache import EventIndex, EventIndexItem


def test_index(tmp_path: Path):
    """Test that we can set, get, and load cached events for an OTU."""
    path = Path(tmp_path / "index")

    index = EventIndex(path)

    otu_id = uuid.uuid4()

    index.set(otu_id, [1, 5, 20, 99, 110], 112)

    assert (
        index.get(otu_id)
        == EventIndex(path).get(otu_id)
        == EventIndexItem(
            112,
            [1, 5, 20, 99, 110],
            otu_id,
        )
    )


def test_index_not_found(tmp_path: Path):
    """Test that we get ``None`` when an OTU ID is not found."""
    path = Path(tmp_path / "index")

    assert EventIndex(path).get(uuid.uuid4()) is None
