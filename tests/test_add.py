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

    subprocess.run([*command, *accessions])


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
    @pytest.mark.parametrize("taxid", [438782])
    def test_empty_success(
        self, taxid, precached_repo: Repo, snapshot: SnapshotAssertion
    ):
        run_add_otu_command(taxid=taxid, path=precached_repo.path, autofill=False)

        for otu in precached_repo.iter_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            assert not otu.accession_set

    @pytest.mark.ncbi
    @pytest.mark.parametrize("taxid", [438782])
    def test_autofill_success(
        self, taxid, precached_repo: Repo, snapshot: SnapshotAssertion
    ):
        run_add_otu_command(taxid=taxid, path=precached_repo.path, autofill=True)

        for otu in precached_repo.iter_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            assert otu.accession_set


class TestAddSequences:
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
    def test_success(
        self, taxid, accessions, precached_repo, snapshot: SnapshotAssertion
    ):
        run_add_sequences_command(taxid, accessions, precached_repo.path)

        for otu in precached_repo.iter_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            for isolate in otu.isolates:
                assert isolate.dict() == snapshot(exclude=props("id", "sequences"))

                sequence_dict = {
                    sequence.accession: sequence for sequence in isolate.sequences
                }

                for accession in sorted(isolate.accession_set):
                    assert sequence_dict[accession].dict() == snapshot(
                        exclude=props("id")
                    )


@pytest.mark.ncbi
class TestUpdateOTU:
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
    def test_success_no_exclusions(
        self, taxid: int, accessions: list[str], precached_repo: Repo
    ):
        run_add_otu_command(taxid, path=precached_repo.path, autofill=False)
        run_add_sequences_command(taxid, accessions, path=precached_repo.path)

        otu = list(precached_repo.iter_otus())[0]

        pre_accession_set = otu.accession_set

        assert set(accessions) == pre_accession_set

        run_update_otu_command(taxid, path=precached_repo.path)

        otu = list(precached_repo.iter_otus())[0]

        assert set(otu.accession_set) != pre_accession_set

        assert set(otu.accession_set) > pre_accession_set
