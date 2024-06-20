import pytest

from virtool_cli.ref.snapshot.index import Snapshotter


@pytest.fixture()
def snapshotter(scratch_repo):
    return Snapshotter(scratch_repo.path / ".cache/snapshot")


class TestSnapshotIndex:
    def test_otu_ids(self, scratch_repo, snapshotter):
        true_otu_ids = [otu.id for otu in scratch_repo.get_all_otus(ignore_cache=True)]

        assert snapshotter.otu_ids

        assert len(true_otu_ids) == len(snapshotter.otu_ids)

        for otu_id in true_otu_ids:
            assert otu_id in snapshotter.id_to_taxid

    def test_taxids(self, scratch_repo, snapshotter):
        true_otu_taxids = [
            otu.taxid for otu in scratch_repo.get_all_otus(ignore_cache=True)
        ]

        assert snapshotter.taxids

        assert len(true_otu_taxids) == len(snapshotter.taxids)

        for taxid in true_otu_taxids:
            assert taxid in snapshotter.index_by_taxid

    def test_names(self, scratch_repo, snapshotter):
        true_otu_names = [
            otu.name for otu in scratch_repo.get_all_otus(ignore_cache=True)
        ]

        assert len(true_otu_names) == len(snapshotter.index_by_name)

        for name in true_otu_names:
            assert name in snapshotter.index_by_name

    def test_accessions(self, scratch_repo, snapshotter):
        true_accessions = set()
        for otu in scratch_repo.get_all_otus(ignore_cache=True):
            true_accessions.update(otu.accessions)

        assert true_accessions

        assert true_accessions == snapshotter.accessions


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
    def test_load_otu_by_taxid(
        self,
        taxid: int,
        accessions: list[str],
        scratch_repo,
        snapshotter,
    ):
        scratch_repo.snapshot()

        rehydrated_otu = scratch_repo.get_otu_by_taxid(taxid)

        assert rehydrated_otu

        snapshot_otu = snapshotter.load_by_taxid(taxid)

        assert snapshot_otu

        assert rehydrated_otu.accessions == snapshot_otu.accessions
