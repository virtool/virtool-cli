import shutil

import pytest
from syrupy import SnapshotAssertion

from virtool_cli.ref.snapshot.model import (
    OTUSnapshotOTU,
    OTUSnapshotIsolate,
    OTUSnapshotSequence,
)

from virtool_cli.ref.resources import (
    EventSourcedRepoOTU,
    EventSourcedRepoIsolate,
    EventSourcedRepoSequence,
)


class TestRepoToSnapshotModel:
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
            (
                1169032,
                ["MH200607", "NC_003355", "KJ207375", "MK431779", "AB017504"],
            ),
        ],
    )
    def test_sequence_conversion(
        self,
        taxid: int,
        accessions: list[str],
        scratch_repo,
        snapshot: SnapshotAssertion,
    ):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for accession in accessions:
            original_sequence = otu.get_sequence_by_accession(accession)

            assert type(original_sequence) is EventSourcedRepoSequence

            converted_model = OTUSnapshotSequence(**original_sequence.dict())

            assert converted_model.model_dump() == snapshot

    @pytest.mark.parametrize("taxid", [438782, 1441799, 430059])
    def test_isolate_conversion(self, taxid, scratch_repo, snapshot: SnapshotAssertion):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for isolate in otu.isolates:
            assert type(isolate) is EventSourcedRepoIsolate

            converted_model = OTUSnapshotIsolate(**isolate.dict())

            assert converted_model.model_dump() == snapshot

    @pytest.mark.parametrize("taxid", [438782, 1441799, 430059])
    def test_otu_conversion(self, taxid, scratch_repo, snapshot: SnapshotAssertion):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        assert type(otu) is EventSourcedRepoOTU

        converted_model = OTUSnapshotOTU(**otu.dict())

        assert converted_model.model_dump(by_alias=True) == snapshot
