from pathlib import Path
import shutil
import subprocess

import pytest

from virtool_cli.ref.repo import EventSourcedRepo
from virtool_cli.ref.snapshot import SnapshotCache


@pytest.fixture()
def cache_example_path(test_files_path):
    return test_files_path / "cache_test"


@pytest.fixture
def precached_repo(tmp_path, test_files_path):
    empty_repo_path = tmp_path / "empty_repo"

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

    shutil.copytree(
        (test_files_path / "cache_test"), (cache_path / "ncbi"), dirs_exist_ok=True
    )
    (cache_path / "repo").mkdir(exist_ok=True)

    yield EventSourcedRepo(empty_repo_path)

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

    subprocess.run(command + accessions, check=False)


class TestSnapshotCache:
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
    def test_snapshot_fidelity(self, taxid, accessions, precached_repo):
        snapshotter = SnapshotCache(precached_repo)

        run_add_otu_command(taxid, path=precached_repo.path, autofill=False)

        run_add_sequences_command(taxid, accessions, path=precached_repo.path)

        subprocess.run(["virtool", "ref", "snapshot"])

        rehydrated_otu = precached_repo.get_otu_by_taxid(taxid)

        snapshot_otu = snapshotter.load_otu(taxid)

        assert rehydrated_otu == snapshot_otu
