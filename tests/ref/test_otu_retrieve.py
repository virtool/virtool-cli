from uuid import UUID

import pytest

from virtool_cli.ref.resources import (
    EventSourcedRepoOTU,
    EventSourcedRepoIsolate,
    EventSourcedRepoSequence,
)


class TestSequence:
    @pytest.mark.parametrize(
        "taxid,accessions",
        [
            (
                345184,
                ["DQ178614", "DQ178613", "DQ178610", "DQ178611"],
            ),
        ],
    )
    def test_equivalence(self, taxid, accessions, scratch_repo):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for accession in accessions:
            sequence = otu.get_sequence_by_accession(accession)

            assert type(sequence) is EventSourcedRepoSequence

            sequence_copy = EventSourcedRepoSequence(**sequence.dict())

            assert sequence == sequence_copy


class TestIsolate:
    @pytest.mark.parametrize("taxid", [345184])
    def test_equivalence(self, taxid, scratch_repo):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for isolate in otu.isolates:
            assert type(isolate) is EventSourcedRepoIsolate

            isolate_copy = EventSourcedRepoIsolate(
                uuid=isolate.id,
                name=isolate.name,
                legacy_id=isolate.legacy_id,
                sequences=isolate.sequences,
            )

            assert isolate == isolate_copy


class TestOTU:
    # @pytest.mark.skip()
    @pytest.mark.parametrize("taxid", [345184])
    def test_equivalence(self, taxid, scratch_repo):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        assert type(otu) is EventSourcedRepoOTU

        otu_copy = EventSourcedRepoOTU(
            uuid=otu.id,
            taxid=otu.taxid,
            name=otu.name,
            acronym=otu.acronym,
            molecule=otu.molecule,
            legacy_id=otu.legacy_id,
            schema=otu.schema,
            excluded_accessions=otu.excluded_accessions,
            isolates=otu.isolates,
        )

        assert otu == otu_copy
