import shutil
import subprocess
from pathlib import Path

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from virtool_cli.ref.repo import EventSourcedRepo as Repo


@pytest.fixture()
def cache_example_path(test_files_path):
    return test_files_path / "cache_test"


@pytest.fixture
def empty_repo_path(tmp_path):
    return tmp_path / "empty_repo"


@pytest.fixture
def precached_repo(empty_repo_path, cache_example_path):
    cache_path = empty_repo_path / ".cache"

    subprocess.run(
        [
            "virtool",
            "ref",
            "init",
            "--data-type",
            "genome",
            "--name",
            "empty repo",
            "--organism",
            "virus",
            "--path",
            str(empty_repo_path),
        ],
        check=False,
    )

    shutil.copytree(cache_example_path, (cache_path / "ncbi"), dirs_exist_ok=True)
    (cache_path / "repo").mkdir(exist_ok=True)

    yield Repo(empty_repo_path)

    shutil.rmtree(empty_repo_path)


def run_add_otu_command(taxid: int, path: Path, autofill: bool):
    otu_command = [
        "virtool",
        "ref",
        "otu",
        "create",
        str(taxid),
        "--path",
        str(path),
    ]

    if autofill:
        otu_command.append("--autofill")

    subprocess.run(otu_command, check=False)


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

    subprocess.run(command + accessions)


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


class TestAddOTU:
    @pytest.mark.parametrize("taxid", [345184])
    def test_empty_success(
        self, taxid, precached_repo: Repo, snapshot: SnapshotAssertion
    ):
        run_add_otu_command(taxid=taxid, path=precached_repo.path, autofill=False)

        for otu in precached_repo.get_all_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            assert not otu.accessions

    @pytest.mark.ncbi
    @pytest.mark.parametrize("taxid", [345184])
    def test_autofill_success(
        self, taxid, precached_repo: Repo, snapshot: SnapshotAssertion
    ):
        run_add_otu_command(taxid=taxid, path=precached_repo.path, autofill=True)

        for otu in precached_repo.get_all_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            assert otu.accessions


class TestAddSequences:
    @pytest.mark.parametrize(
        "taxid,accessions",
        [
            (
                345184,
                ["DQ178614", "DQ178613", "DQ178610", "DQ178611"],
            ),
        ],
    )
    def test_success(
        self, taxid, accessions, precached_repo, snapshot: SnapshotAssertion
    ):
        run_add_sequences_command(taxid, accessions, precached_repo.path)

        for otu in precached_repo.get_all_otus():
            assert otu.accessions == set(accessions)

            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            for isolate in otu.isolates:
                assert isolate.dict() == snapshot(exclude=props("id", "sequences"))

                for accession in sorted(isolate.accessions):
                    assert isolate.get_sequence_by_accession(
                        accession
                    ).dict() == snapshot(exclude=props("id"))


@pytest.mark.ncbi
class TestUpdateOTU:
    @pytest.mark.parametrize(
        "taxid,accessions",
        [
            (
                345184,
                ["DQ178614", "DQ178613", "DQ178610", "DQ178611"],
            ),
        ],
    )
    def test_success_no_exclusions(
        self, taxid: int, accessions: list[str], precached_repo: Repo
    ):
        run_add_otu_command(taxid, path=precached_repo.path, autofill=False)
        run_add_sequences_command(taxid, accessions, path=precached_repo.path)

        otu = precached_repo.get_all_otus()[0]

        pre_accession_set = otu.accessions

        assert set(accessions) == pre_accession_set

        run_update_otu_command(taxid, path=precached_repo.path)

        otu = precached_repo.get_all_otus()[0]

        assert set(otu.accessions) != pre_accession_set

        assert set(otu.accessions) > pre_accession_set
