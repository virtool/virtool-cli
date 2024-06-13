import shutil
from typing import List, Dict
from uuid import UUID

import orjson
import pytest
from pydantic import TypeAdapter, ValidationError
from syrupy import SnapshotAssertion

from virtool_cli.ref.snapshot.otu import OTUSnapshot
from virtool_cli.ref.snapshot.model import (
    OTUSnapshotOTU,
    OTUSnapshotIsolate,
    OTUSnapshotToCIsolate,
)


@pytest.fixture()
def empty_otu_snapshot_path(tmp_path):
    return tmp_path / "empty_otu_snapshot"


@pytest.fixture()
def snapshotted_scratch(scratch_repo):
    scratch_repo.snapshot()

    yield scratch_repo

    shutil.rmtree(scratch_repo.cache_path / "snapshot")


def is_accession_in_toc(accession: str, toc: dict) -> bool:
    for key in toc:
        isolate = toc[key]
        if accession in isolate["accessions"]:
            return True

    return False


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

    def test_otu_data(self, taxid: int, snapshotted_scratch):
        otu_id = snapshotted_scratch.index_by_taxid[taxid]
        otu_snapshot_path = snapshotted_scratch.cache_path / f"snapshot/{otu_id}"

        with open((otu_snapshot_path / "otu.json"), "rb") as f:
            otu_model = OTUSnapshotOTU.model_validate_json(f.read())

        assert otu_model.id == otu_id

        assert otu_model.taxid == taxid

    def test_toc(self, taxid: int, snapshotted_scratch, snapshot):
        otu_id = snapshotted_scratch.index_by_taxid[taxid]
        otu_snapshot_path = snapshotted_scratch.cache_path / f"snapshot/{otu_id}"

        with open((otu_snapshot_path / "toc.json"), "rb") as f:
            toc_dict = orjson.loads(f.read())

        toc_adapter = TypeAdapter(Dict[str, OTUSnapshotToCIsolate])

        try:
            toc_adapter.validate_python(toc_dict)
        except ValidationError as e:
            print(e)

        rehydrated_otu = snapshotted_scratch.get_otu(otu_id)

        for accession in rehydrated_otu.accessions:
            assert is_accession_in_toc(accession, toc_dict)
