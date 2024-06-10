from uuid import UUID

import pytest

from virtool_cli.ref.snapshot.otu import OTUSnapshot


@pytest.fixture()
def empty_otu_snapshot_path(tmp_path):
    return tmp_path / "empty_otu_snapshot"


class TestOTUSnapshot:
    @pytest.mark.parametrize("taxid", [438782, 1441799, 430059])
    def test_cache(self, taxid: int, scratch_repo, empty_otu_snapshot_path):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        otu_snapshotter = OTUSnapshot(empty_otu_snapshot_path)

        otu_snapshotter.cache(otu, options=None)

        assert (otu_snapshotter.path / "otu.json").exists()

        assert (otu_snapshotter.path / "toc.json").exists()

        data_dir = otu_snapshotter.path / "data"

        assert data_dir.exists()

        data_set = set(UUID(data_path.stem) for data_path in data_dir.glob("*.json"))

        isolate_ids = set(isolate.id for isolate in otu.isolates)

        assert isolate_ids.issubset(data_set)
