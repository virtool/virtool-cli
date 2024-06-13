from pathlib import Path
import shutil
import subprocess

import pytest


class TestSnapshotIndex:
    def test_otu_ids(self, scratch_repo):
        true_otu_ids = [otu.id for otu in scratch_repo.get_all_otus(ignore_cache=True)]

        assert scratch_repo._snapshotter.otu_ids

        assert len(true_otu_ids) == len(scratch_repo._snapshotter.otu_ids)

        for otu_id in true_otu_ids:
            assert otu_id in scratch_repo._snapshotter.index_by_id

    def test_taxids(self, scratch_repo):
        true_otu_taxids = [
            otu.taxid for otu in scratch_repo.get_all_otus(ignore_cache=True)
        ]

        assert scratch_repo._snapshotter.taxids

        assert len(true_otu_taxids) == len(scratch_repo._snapshotter.taxids)

        for taxid in true_otu_taxids:
            assert taxid in scratch_repo._snapshotter.index_by_taxid

    def test_names(self, scratch_repo):
        true_otu_names = [
            otu.name for otu in scratch_repo.get_all_otus(ignore_cache=True)
        ]

        assert scratch_repo._snapshotter.taxids

        assert len(true_otu_names) == len(scratch_repo._snapshotter.index_by_name)

        for name in true_otu_names:
            assert name in scratch_repo._snapshotter.index_by_name


class TestSnapshotIndexCaching:
    @pytest.mark.parametrize(
        "taxid,accessions",
        [
            (
                438782,
                [
                    "NC_010314",
                    "NC_010315",
                    "NC_010316",
                    "NC_010317",
                    "NC_010318",
                    "NC_010319",
                ],
            ),
        ],
    )
    def test_load_otu_by_taxid(self, taxid: int, accessions: list[str], scratch_repo):
        scratch_repo.snapshot()

        rehydrated_otu = scratch_repo.get_otu_by_taxid(taxid)

        assert rehydrated_otu

        snapshot_otu = scratch_repo._snapshotter.load_otu_by_taxid(taxid)

        assert snapshot_otu

        assert rehydrated_otu.accessions == snapshot_otu.accessions
