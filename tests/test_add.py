import shutil
import subprocess
from pathlib import Path

import orjson
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


def index_repo_data(repo: Repo):
    otus = {}

    for otu in repo.iter_otus():
        otus[otu.taxid] = otu.dict()

        isolates = {}
        for isolate in otu.isolates:
            isolates[(isolate.name.type, isolate.name.value)] = isolate.dict()

            sequences = {}
            for sequence in isolate.sequences:
                sequences[sequence.accession] = sequence.dict()

            isolates[(isolate.name.type, isolate.name.value)]["sequences"] = sequences

        otus[otu.taxid]["isolates"] = isolates

    return otus


class TestAddOTU:
    @staticmethod
    def run_command(taxid: int, path: Path, autofill: bool):
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

    @pytest.mark.parametrize("taxid", [438782])
    def test_empty_success(
        self, taxid, precached_repo: Repo, snapshot: SnapshotAssertion
    ):
        self.run_command(taxid=taxid, path=precached_repo.path, autofill=False)

        for otu in precached_repo.iter_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            assert not otu.accessions

    @pytest.mark.ncbi
    @pytest.mark.parametrize("taxid", [438782])
    def test_autofill_success(
        self, taxid, precached_repo: Repo, snapshot: SnapshotAssertion
    ):
        self.run_command(taxid=taxid, path=precached_repo.path, autofill=True)

        for otu in precached_repo.iter_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            assert otu.accessions


class TestAddSequences:
    @staticmethod
    def run_add_sequences(taxid: int, accessions: list[str], path: Path):
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
        self.run_add_sequences(taxid, accessions, precached_repo.path)

        for otu in precached_repo.iter_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            for isolate in otu.isolates:
                assert isolate.dict() == snapshot(exclude=props("id", "sequences"))

                sequence_dict = {
                    sequence.accession: sequence for sequence in isolate.sequences
                }

                for accession in sorted(isolate.accessions):
                    assert sequence_dict[accession].dict() == snapshot(
                        exclude=props("id")
                    )
