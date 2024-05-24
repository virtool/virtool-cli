import pytest
import shutil
import subprocess

from syrupy import SnapshotAssertion
from syrupy.filters import props


@pytest.fixture
def precached_repo(empty_repo, test_files_path):
    shutil.copytree(
        (test_files_path / "cache_test"),
        (empty_repo.path / ".cache/ncbi"),
        dirs_exist_ok=True,
    )

    yield empty_repo

    shutil.rmtree(empty_repo.path)


@pytest.fixture
def ref_table_of_contents() -> dict[int, list[str]]:
    return {
        96892: ["NC_005954"],
        683179: ["NC_013464"],
        438782: [
            "NC_010314",
            "NC_010318",
            "NC_010316",
            "NC_010319",
            "NC_010315",
            "NC_010317",
            "EF546804.1",
            "EF546802.1",
            "EF546803.1",
            "EF546806.1",
            "EF546807.1",
            "EF546805.1",
        ],
        1169032: ["MH200607.1", "NC_003355", "KJ207375.1", "MK431779.1", "AB017504.1"],
        270478: ["NC_011560"],
        1468454: ["NC_023641"],
        1278205: ["NC_020160"],
        1505530: ["NC_024301"],
        48201: ["NC_043133.1"],
        1468172: ["NC_023892"],
        946046: ["NC_014791"],
        1441799: ["NC_023881"],
        132477: ["NC_013006"],
        1414644: ["NC_022745"],
        56879: ["NC_001793"],
        223262: ["NC_004630", "NC_004625"],
        2060511: ["NC_036587"],
        1324128: ["NC_021786"],
        345184: ["DQ178614", "DQ178613", "DQ178610", "DQ178611"],
        430059: ["NC_043170.1"],
        2170135: ["NC_030236.1"],
        1561150: ["NC_026472"],
        1198450: ["NC_038796.1", "NC_038797"],
        662596: ["NC_010239"],
    }


@pytest.mark.ncbi()
def test_add_ref_mini(
    ref_table_of_contents, precached_repo, snapshot: SnapshotAssertion
):
    for taxid in ref_table_of_contents:
        subprocess.run(
            [
                "virtool",
                "ref",
                "otu",
                "create",
                "--path",
                str(precached_repo.path),
                str(taxid),
            ],
            check=False,
        )

        subprocess.run(
            [
                "virtool",
                "ref",
                "sequences",
                "add",
                "--path",
                str(precached_repo.path),
                "--taxid",
                str(taxid),
            ]
            + ref_table_of_contents[taxid],
            check=False,
        )

    for taxid in ref_table_of_contents:
        otu = precached_repo.get_otu_by_taxid(taxid)
        assert otu.dict() == snapshot(exclude=props("id", "isolates"))
        assert otu.accessions == snapshot

        assert otu.accessions == set(
            accession.split(".")[0] for accession in ref_table_of_contents[taxid]
        )
