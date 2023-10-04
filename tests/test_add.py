import pytest
from pathlib import Path
import json
import shutil
import subprocess

from virtool_cli.utils.reference import get_isolate_paths

from paths import TEST_FILES_PATH

BASE_PATH = TEST_FILES_PATH / "src_test"
TEST_ACCLOG_PATH = TEST_FILES_PATH / "catalog"
TEST_ACCESSION = "NC_010319"


@pytest.fixture()
def work_path(tmp_path):
    test_src_path = tmp_path / "src"
    shutil.copytree(BASE_PATH, test_src_path)
    return test_src_path


class TestAddAccession:
    @staticmethod
    def run_command(accession: str, src_path: Path, catalog_path: Path):
        subprocess.call(
            [
                "virtool",
                "ref",
                "add",
                "accession",
                "-acc",
                accession,
                "-src",
                str(src_path),
                "-cat",
                str(catalog_path),
            ]
        )

    def run_add_accession(
        self, accession: str, otu_dirname: str, src_path: Path, catalog_path: Path
    ):
        """
        Add into an existing isolate directory
        """
        isolate_path = src_path / otu_dirname

        pre_sequence_paths = set(isolate_path.glob("*.json"))

        self.run_command(accession, src_path, catalog_path)

        post_sequence_paths = set(isolate_path.glob("*.json"))

        return pre_sequence_paths, post_sequence_paths

    def run_add_isolate(
        self, accession: str, otu_subpath: str, src_path: Path, catalog_path: Path
    ):
        """
        Add into a new isolate directory
        """
        otu_path = src_path / otu_subpath

        pre_isolate_paths = set(get_isolate_paths(otu_path))

        self.run_command(accession, src_path, catalog_path)

        post_isolate_paths = set(get_isolate_paths(otu_path))

        return pre_isolate_paths, post_isolate_paths

    @pytest.mark.parametrize(
        "accession, isolate_subpath",
        [("DQ178612", "cabbage_leaf_curl_jamaica_virus--d226290f/d293d531")],
    )
    def test_add_accession_success(self, accession, isolate_subpath, work_path):
        """
        Check that virtool ref add accession does the job when the isolate exists
        """
        pre_sequence_paths, post_sequence_paths = self.run_add_accession(
            accession,
            isolate_subpath,
            src_path=work_path,
            catalog_path=TEST_ACCLOG_PATH,
        )

        new_sequences = post_sequence_paths.difference(pre_sequence_paths)

        assert new_sequences

    @pytest.mark.parametrize(
        "accession, isolate_subpath",
        [("NC_010319", "abaca_bunchy_top_virus--c93ec9a9/4e8amg20")],
    )
    def test_add_accession_fail(self, accession, isolate_subpath, work_path):
        """
        Check that virtool ref add accession does not add sequences that already exist
        """
        pre_sequence_paths, post_sequence_paths = self.run_add_accession(
            accession,
            isolate_subpath,
            src_path=work_path,
            catalog_path=TEST_ACCLOG_PATH,
        )

        new_sequences = post_sequence_paths.difference(pre_sequence_paths)

        assert not new_sequences

    @pytest.mark.parametrize(
        "accession, otu_dirname",
        [("KT390494", "nanovirus-like_particle--ae0f2a35")],
    )
    def test_add_isolate_success(self, accession, otu_dirname, work_path):
        """
        Check that virtool ref add accession does the job when the isolate does not exist
        """
        pre_isolate_paths, post_isolate_paths = self.run_add_isolate(
            accession, otu_dirname, src_path=work_path, catalog_path=TEST_ACCLOG_PATH
        )

        new_isolates = post_isolate_paths.difference(pre_isolate_paths)

        print(pre_isolate_paths)

        print(post_isolate_paths)

        assert new_isolates

        assert (new_isolates.pop() / "isolate.json").exists()


class TestAddOTU:
    @staticmethod
    def run_command(taxon_id: int, src_path: Path, catalog_path: Path):
        subprocess.call(
            [
                "virtool",
                "ref",
                "add",
                "otu",
                "-taxid",
                str(taxon_id),
                "-src",
                str(src_path),
                "-cat",
                str(catalog_path),
            ]
        )

    def run_add_otu(self, taxon_id: int, src_path: Path, catalog_path: Path):
        """
        Attempt to add a new OTU
        """

        pre_otu_paths = set(src_path.glob("*--*"))

        self.run_command(taxon_id, src_path, catalog_path)

        post_otu_paths = set(src_path.glob("*--*"))

        return pre_otu_paths, post_otu_paths

    @pytest.mark.parametrize("taxon_id", [908125])
    def test_add_otu_success(self, taxon_id, work_path):
        """ """
        pre_otu_paths = set(work_path.glob("*--*"))

        self.run_add_otu(
            taxon_id=taxon_id, src_path=work_path, catalog_path=TEST_ACCLOG_PATH
        )

        post_otu_paths = set(work_path.glob("*--*"))

        new_otus = post_otu_paths.difference(pre_otu_paths)

        assert new_otus

        otu_path = new_otus.pop()

        assert (otu_path / "otu.json").exists()

    @pytest.mark.parametrize("taxon_id", [345184])
    def test_add_otu_conflict(self, taxon_id, work_path):
        """Don't add taxon IDs that already exist"""
        pre_otu_paths = set(work_path.glob("*--*"))

        self.run_add_otu(
            taxon_id=taxon_id, src_path=work_path, catalog_path=TEST_ACCLOG_PATH
        )

        post_otu_paths = set(work_path.glob("*--*"))

        new_otus = post_otu_paths.difference(pre_otu_paths)

        assert not new_otus
