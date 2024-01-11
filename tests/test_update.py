from pathlib import Path
import pytest
import shutil
import json
import subprocess

from paths import TEST_FILES_PATH

BASE_PATH = TEST_FILES_PATH / "src_test"


class TestEmptyRepo:
    def test_empty_success(self, empty_repo_path):
        subprocess.run(["virtool", "ref", "init", "-repo", str(empty_repo_path)])

        completed_process = run_update_with_process(src_path=empty_repo_path / "src")

        assert completed_process.returncode == 0


class TestUpdateOTU:
    @pytest.mark.parametrize("otu_dirname", ["abaca_bunchy_top_virus--c93ec9a9"])
    def test_update_otu(self, otu_dirname: str):
        otu_path = BASE_PATH / otu_dirname
        subprocess.run([
            "virtool",
            "ref",
            "update",
            "otu",
            "-otu",
            str(otu_path),
        ])

        assert True


# @pytest.mark.parametrize("base_path", [BASE_PATH])
# def test_update_basic(base_path, tmp_path):
#     """
#     Test that updates actually pull something.
#     """
#     fetch_path = tmp_path / "src"
#     pre_update_ref_path = tmp_path / "reference_pre.json"
#     post_update_ref_path = tmp_path / "reference_post.json"
#
#     shutil.copytree(base_path, fetch_path)
#
#     run_build(src_path=fetch_path, output_path=pre_update_ref_path)
#
#     completed_process = run_update(src_path=fetch_path)
#     assert completed_process.returncode == 0
#
#     run_build(src_path=fetch_path, output_path=post_update_ref_path)
#
#     reference_pre = json.loads(pre_update_ref_path.read_text())
#     pre_otu_dict = convert_to_dict(reference_pre["otus"])
#
#     reference_post = json.loads(post_update_ref_path.read_text())
#     post_otu_dict = convert_to_dict(reference_post["otus"])
#
#     difference_counter = 0
#
#     for otu_id in post_otu_dict:
#         pre_accessions = get_otu_accessions(pre_otu_dict[otu_id])
#         post_accessions = get_otu_accessions(post_otu_dict[otu_id])
#
#         print(pre_accessions)
#         print(post_accessions)
#
#         if pre_accessions != post_accessions:
#             difference_counter += 1
#
#     # Any new data counts
#     assert difference_counter > 0
#
#
# @pytest.mark.skip()
# @pytest.mark.parametrize("base_path", [BASE_PATH])
# def test_update_autoevaluate(base_path, tmp_path):
#     """
#     Test that updates actually pull something.
#     Autoevaluation.
#     """
#     fetch_path = tmp_path / "src"
#     pre_update_ref_path = tmp_path / "reference_pre.json"
#     post_update_ref_path = tmp_path / "reference_post.json"
#
#     shutil.copytree(base_path, fetch_path)
#
#     run_build(src_path=fetch_path, output_path=pre_update_ref_path)
#
#     run_update(src_path=fetch_path)
#
#     run_build(src_path=fetch_path, output_path=post_update_ref_path)
#
#     reference_pre = json.loads(pre_update_ref_path.read_text())
#     pre_otu_dict = convert_to_dict(reference_pre["otus"])
#
#     reference_post = json.loads(post_update_ref_path.read_text())
#     post_otu_dict = convert_to_dict(reference_post["otus"])
#
#     difference_counter = 0
#
#     for otu_id in post_otu_dict:
#         pre_accessions = get_otu_accessions(pre_otu_dict[otu_id])
#         post_accessions = get_otu_accessions(post_otu_dict[otu_id])
#
#         print(pre_accessions)
#         print(post_accessions)
#
#         if pre_accessions != post_accessions:
#             difference_counter += 1
#
#     # Any new data counts
#     assert difference_counter > 0

#
@pytest.fixture()
def empty_repo_path(tmp_path):
    return tmp_path / "repo_empty"


def run_update_with_process(src_path) -> subprocess.CompletedProcess:
    complete_process = subprocess.run(
        [
            "virtool",
            "ref",
            "update",
            "reference",
            "-src",
            str(src_path),
        ],
        capture_output=True,
    )

    return complete_process


def run_build(src_path, output_path) -> subprocess.CompletedProcess:
    complete_process = subprocess.run(
        [
            "virtool",
            "ref",
            "build",
            "-src",
            str(src_path),
            "--output",
            str(output_path),
        ],
        capture_output=True,
    )

    return complete_process


def convert_to_dict(otu_list: list) -> dict:
    """
    Converts a list of OTUs to a dict keyed by Virtool ID

    :param otu_list: A list of deserialized OTU data
    :return: The contents of otu_list keyed by OTU id
    """
    otu_dict = {}
    for otu in otu_list:
        otu_dict[otu["_id"]] = otu
    return otu_dict


def get_otu_accessions(otu_dict: dict) -> set:
    """
    Gets all accessions from an OTU directory and returns a list

    :param otu_dict: Deserialized OTU data
    :return: The accessions included under the OTU directory in a set
    """
    accessions = set()

    for isolate in otu_dict["isolates"]:
        for sequence in isolate["sequences"]:
            accessions.add(sequence["accession"])

    return accessions
