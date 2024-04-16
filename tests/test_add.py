import shutil
import subprocess
from pathlib import Path

import orjson
import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from virtool_cli.ref.repo import EventSourcedRepo as Repo
from virtool_cli.utils.reference import get_sequence_paths


@pytest.fixture
def empty_repo_path(tmp_path):
    yield tmp_path / "empty_repo"


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


def get_all_sequence_paths(otu_path: Path) -> set[Path]:
    return {
        sequence_path
        for isolate_path in otu_path.iterdir()
        if isolate_path.is_dir()
        for sequence_path in get_sequence_paths(isolate_path)
    }


def run_build_reference(repo: Path, output: Path):
    subprocess.run(
        [
            "virtool",
            "ref",
            "build",
            "--path",
            str(repo),
            "--output-path",
            str(output),
        ],
        check=False,
    )


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
    def run_add_otu(taxid: int, accessions: list[str], path: Path):
        command = [
            "virtool",
            "ref",
            "otu",
            "init",
            "--path",
            str(path),
            "--taxid",
            str(taxid),
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
        self.run_add_otu(taxid, accessions, precached_repo.path)

        for otu in precached_repo.iter_otus():
            assert otu.dict() == snapshot(exclude=props("id", "isolates"))

            assert sorted(otu.get_indexed_isolates()) == snapshot

            assert sorted(otu.get_accessions()) == snapshot


@pytest.mark.ncbi
class TestAddOTUAutofill:
    @staticmethod
    def run_add_otu(taxid: int, path: Path):
        subprocess.run(
            ["virtool", "ref", "otu", "create", str(taxid), "--path", str(path)],
            check=False,
        )

    @staticmethod
    def run_build_reference(repo: Path, output: Path):
        subprocess.run(
            [
                "virtool",
                "ref",
                "build",
                "--path",
                str(repo),
                "--output-path",
                str(output),
            ],
            check=False,
        )

    @staticmethod
    def read_build(path) -> dict:
        with open(path, "rb") as f:
            return orjson.loads(f.read())

    @pytest.mark.parametrize("taxid", [438782])
    def test_success(self, taxid, precached_repo: Repo):
        print(list((precached_repo.path / ".cache").iterdir()))
        output_path = precached_repo.path / "reference.json"

        self.run_build_reference(precached_repo.path, precached_repo.path)
        pre_repo = self.read_build(output_path)

        assert pre_repo["otus"] == []

        self.run_add_otu(taxid=taxid, path=precached_repo.path)

        self.run_build_reference(precached_repo.path, precached_repo.path)

        assert output_path.exists()
        post_repo = self.read_build(output_path)

        assert post_repo["otus"]

        for key in ["created_at", "data_type", "name"]:
            assert pre_repo[key] == post_repo[key]

        assert post_repo["otus"] != pre_repo["otus"]
