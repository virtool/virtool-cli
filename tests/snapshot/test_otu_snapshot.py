from uuid import UUID

import orjson
import pytest

from virtool_cli.ref.snapshot.otu import OTUSnapshot


@pytest.fixture()
def empty_otu_snapshot_path(tmp_path):
    return tmp_path / "empty_otu_snapshot"


@pytest.mark.parametrize("taxid", [438782, 1441799, 430059])
class TestOTUSnapshot:
    def test_snapshot_direct(self, taxid: int, scratch_repo, empty_otu_snapshot_path):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        otu_snapshotter = OTUSnapshot(empty_otu_snapshot_path)

        otu_snapshotter.cache(otu, options=None)

        assert otu_snapshotter.at_event is None

        assert (otu_snapshotter.path / "otu.json").exists()

        assert (otu_snapshotter.path / "toc.json").exists()

        assert (otu_snapshotter.path / "metadata.json").exists()

        data_dir = otu_snapshotter.path / "data"

        assert data_dir.exists()

        data_set = set(UUID(data_path.stem) for data_path in data_dir.glob("*.json"))

        isolate_ids = set(isolate.id for isolate in otu.isolates)

        assert isolate_ids.issubset(data_set)

    def test_snapshot_from_index(self, taxid: int, scratch_repo):
        scratch_repo.snapshot()

        otu_id = scratch_repo.index_by_taxid[taxid]

        otu_snapshot_path = scratch_repo.cache_path / f"snapshot/{otu_id}"

        assert otu_snapshot_path.exists()

        assert otu_snapshot_path / "metadata.json"

        with open((otu_snapshot_path / "metadata.json"), "rb") as f:
            metadata_dict = orjson.loads(f.read())

        assert type(metadata_dict["at_event"]) is int
