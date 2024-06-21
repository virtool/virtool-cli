import os
import subprocess
from pathlib import Path

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from virtool_cli.ref.otu import create_otu, update_otu
from virtool_cli.ref.repo import EventSourcedRepo


def run_add_sequences_command(taxid: int, accessions: list[str], path: Path):
    command = [
        "virtool",
        "ref",
        "sequences",
        "add",
        "--taxid",
        str(taxid),
        "--path",
        str(path),
    ]

    subprocess.run(command + accessions, check=False)


def run_update_otu_command(taxid: int, path: Path):
    otu_command = [
        "virtool",
        "ref",
        "otu",
        "update",
        str(taxid),
        "--path",
        str(path),
    ]

    subprocess.run(otu_command, check=False)


class TestCreateOTU:
    def test_empty_success(
        self,
        precached_repo: EventSourcedRepo,
        snapshot: SnapshotAssertion,
    ):
        """Test that an OTU can be created in an empty repository."""
        print(os.listdir(precached_repo.path / "src"))

        otu = create_otu(precached_repo, 345184)

        assert otu.dict() == snapshot(exclude=props("id", "isolates"))

        # Ensure only one OTU is present in the repository, and it matches the return
        # value of the creation function.
        assert list(precached_repo.iter_otus()) == [otu]

    @pytest.mark.ncbi()
    def test_autofill(
        self,
        precached_repo: EventSourcedRepo,
        snapshot: SnapshotAssertion,
    ):
        subprocess.run(
            [
                "virtool",
                "ref",
                "otu",
                "create",
                "345184",
                "--path",
                str(precached_repo.path),
                "--autofill",
            ],
            check=False,
        )

        otus = list(EventSourcedRepo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        otu = otus[0]

        assert otu.dict() == snapshot(exclude=props("id", "isolates"))
        assert otu.accessions


class TestAddSequences:
    def test_success(
        self,
        precached_repo: EventSourcedRepo,
        snapshot: SnapshotAssertion,
    ):
        accessions = ["DQ178614", "DQ178613", "DQ178610", "DQ178611"]
        run_add_sequences_command(345184, accessions, precached_repo.path)

        for otu in precached_repo.get_all_otus():
            assert otu.accessions == set(accessions)

            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            for isolate in otu.isolates:
                assert isolate.dict() == snapshot(exclude=props("id", "sequences"))

                for accession in sorted(isolate.accessions):
                    assert isolate.get_sequence_by_accession(
                        accession,
                    ).dict() == snapshot(exclude=props("id"))


@pytest.mark.ncbi()
class TestUpdateOTU:
    def test_success_no_exclusions(
        self,
        precached_repo: EventSourcedRepo,
        snapshot: SnapshotAssertion,
    ):
        otu = create_otu(precached_repo, 345184)
        update_otu(precached_repo, otu)

        otu = precached_repo.get_otu(otu.id)

        assert [otu.dict() for otu in precached_repo.iter_otus()] == snapshot(
            exclude=props("id", "isolates"),
        )

        assert otu.accessions == {
            "DQ178608",
            "DQ178609",
            "DQ178610",
            "DQ178611",
            "DQ178612",
            "DQ178613",
            "DQ178614",
            "NC_038792",
            "NC_038793",
        }
