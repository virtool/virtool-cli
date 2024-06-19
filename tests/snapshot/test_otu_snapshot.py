import shutil
from uuid import UUID

import orjson
import pytest

from virtool_cli.ref.snapshot.otu import OTUSnapshot


@pytest.fixture()
def empty_otu_snapshot_path(tmp_path):
    return tmp_path / "empty_otu_snapshot"


@pytest.mark.parametrize(
    "otu_id",
    [
        # 438782
        UUID("1f88507e-1a9a-4ec7-88ed-cc253c6200ce"),
        # 1441799
        UUID("c2a22914-12cf-4b66-aa42-e1395e4498de"),
        # 430059
        UUID("352325fa-bff2-44b1-a5ad-64de5c786c74"),
    ],
)
class TestOTUSnapshot:
    def test_snapshot_direct(self, otu_id: UUID, scratch_repo, empty_otu_snapshot_path):
        """Test that OTUSnapshot can build a snapshot directly from a floating OTU."""

        otu = scratch_repo.get_otu(otu_id)

        otu_snapshotter = OTUSnapshot(empty_otu_snapshot_path)

        otu_snapshotter.cache(otu, indent=None)

        assert otu_snapshotter.at_event is None

        for filestem in ("otu", "toc", "metadata"):
            assert (otu_snapshotter.path / f"{filestem}.json").exists()

        data_path = otu_snapshotter.path / "data"

        assert data_path.exists()

        data_set = {UUID(path.stem) for path in data_path.glob("*.json")}

        isolate_ids = {isolate.id for isolate in otu.isolates}

        assert isolate_ids.issubset(data_set)

    def test_snapshot_metadata(self, otu_id: UUID, scratch_repo):
        """Test that EventSourcedRepo's snapshot function produces proper OTU metadata (at_event)."""
        otu_snapshot_path = scratch_repo.cache_path / f"snapshot/{otu_id}"

        assert otu_snapshot_path.exists()

        assert otu_snapshot_path / "metadata.json"

        with open((otu_snapshot_path / "metadata.json"), "rb") as f:
            metadata_dict = orjson.loads(f.read())

        assert type(metadata_dict["at_event"]) is int

    def test_load(self, otu_id: UUID, scratch_repo):
        """Test OTUSnapshot.load()"""

        rehydrated_otu = scratch_repo.get_otu(otu_id)

        assert rehydrated_otu

        otu_snapshotter = OTUSnapshot(
            path=scratch_repo.cache_path / f"snapshot/{otu_id}"
        )

        snapshot_otu = otu_snapshotter.load()

        assert snapshot_otu

        assert rehydrated_otu.accessions == snapshot_otu.accessions

    def test_toc(self, otu_id: UUID, scratch_repo):
        """Test that the table of contents is written correctly."""

        rehydrated_otu = scratch_repo.get_otu(otu_id)

        otu_snapshotter = OTUSnapshot(
            path=scratch_repo.cache_path / f"snapshot/{otu_id}"
        )

        with open((otu_snapshotter.path / "toc.json"), "rb") as f:
            toc_dict = orjson.loads(f.read())

        for isolate in rehydrated_otu.isolates:
            key = str(isolate.name)

            assert key in toc_dict

            for accession in isolate.accessions:
                assert accession in toc_dict[key]["accessions"]
