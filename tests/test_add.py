import subprocess
from pathlib import Path

import pytest

from virtool_cli.ref.build import build_json
from virtool_cli.ref.init import init_reference
from virtool_cli.utils.reference import get_isolate_paths, get_sequence_paths


def get_all_sequence_paths(otu_path: Path) -> set[Path]:
    sequence_paths = set()

    for isolate_path in otu_path.iterdir():
        if isolate_path.is_dir():
            for sequence_path in get_sequence_paths(isolate_path):
                sequence_paths.add(sequence_path)

    return sequence_paths


class TestAddAccession:
    @staticmethod
    def run_add_accession(accession: str, path: Path):
        subprocess.run(
            [
                "virtool",
                "ref",
                "add",
                "accession",
                "--path",
                str(path),
                accession,
            ],
            check=True,
        )

    @pytest.mark.parametrize(
        "accession",
        ["DQ178612", "NC_038793"],
    )
    def test_success(
        self,
        accession: str,
        scratch_path: Path,
    ):
        """Test that an accession is added when the parent isolate already exists."""
        otu_path = scratch_path / "src" / "cabbage_leaf_curl_jamaica_virus--d226290f"

        pre_sequence_paths = get_all_sequence_paths(otu_path)

        self.run_add_accession(accession, scratch_path)

        post_sequence_paths = get_all_sequence_paths(otu_path)

        assert post_sequence_paths - pre_sequence_paths

    @pytest.mark.parametrize(
        "accession, otu_dirname",
        [("KT390494", "nanovirus_like_particle--ae0f2a35")],
    )
    def test_success_new_isolate(self, accession, otu_dirname, scratch_path: Path):
        """Test that an accession is added when the parent isolate does not exist"""
        otu_path = scratch_path / "src" / otu_dirname

        pre_isolate_paths = set(get_isolate_paths(otu_path))

        self.run_add_accession(accession, scratch_path)

        new_isolates = set(get_isolate_paths(otu_path)) - pre_isolate_paths

        assert len(new_isolates) == 1
        assert (new_isolates.pop() / "isolate.json").exists()

    def test_invalid_accession(self, scratch_path: Path):
        otu_path = (
            scratch_path / "src" / "pagoda_yellow_mosaic_associated_virus--dd21fd8f"
        )

        pre_sequence_paths = get_all_sequence_paths(otu_path)

        self.run_add_accession("NC-024301", scratch_path)

        assert get_all_sequence_paths(otu_path) == pre_sequence_paths

    def test_add_accession_fail(
        self,
        scratch_path: Path,
    ):
        """Test that sequences cannot be added that already exist."""
        otu_path = scratch_path / "src" / "abaca_bunchy_top_virus--c93ec9a9"

        pre_sequence_paths = get_all_sequence_paths(otu_path)

        self.run_add_accession("NC_010319", scratch_path)

        post_sequence_paths = get_all_sequence_paths(otu_path)

        assert post_sequence_paths == pre_sequence_paths


class TestAddAccessions:
    @staticmethod
    def run_add_accessions(accessions: list[str], otu_dirname: str, path: Path):
        """Add into an existing isolate directory"""
        otu_path = path / "src" / otu_dirname

        pre_sequence_paths = get_all_sequence_paths(otu_path)

        subprocess.run(
            [
                "virtool",
                "ref",
                "add",
                "accessions",
                "--otu-path",
                str(otu_path),
                *accessions,
            ],
            check=False,
        )

        return pre_sequence_paths, get_all_sequence_paths(otu_path)

    @pytest.mark.parametrize(
        "accessions, otu_dirname",
        [
            (["DQ178612, NC_038793"], "cabbage_leaf_curl_jamaica_virus--d226290f"),
            (["KT390494, KT390496, KT390501"], "nanovirus_like_particle--ae0f2a35"),
        ],
    )
    def test_success(
        self,
        accessions: list[str],
        otu_dirname: str,
        scratch_path: Path,
    ):
        """Check that virtool ref add accessions does the job"""
        pre_sequence_paths, post_sequence_paths = self.run_add_accessions(
            accessions,
            otu_dirname,
            scratch_path,
        )

        assert post_sequence_paths - pre_sequence_paths

    @pytest.mark.parametrize(
        "accessions, otu_dirname",
        [
            (["DQ178611, DQ178610"], "cabbage_leaf_curl_jamaica_virus--d226290f"),
            (["KTX390494, KTX390496, KTX390501"], "nanovirus_like_particle--ae0f2a35"),
        ],
    )
    def test_failure(
        self,
        accessions: list[str],
        otu_dirname: str,
        scratch_path: Path,
    ):
        """Check that virtool ref add accessions does the job"""
        pre_sequence_paths, post_sequence_paths = self.run_add_accessions(
            accessions,
            otu_dirname,
            scratch_path,
        )

        assert post_sequence_paths == pre_sequence_paths


class TestAddOTU:
    @staticmethod
    def run_add_otu(path: Path, taxid: int):
        subprocess.run(
            [
                "virtool",
                "ref",
                "add",
                "otu",
                "--path",
                str(path),
                str(taxid),
            ],
            check=False,
        )

    def test_success(self, scratch_path: Path):
        """ """
        pre_otu_paths = set((scratch_path / "src").glob("*--*"))

        self.run_add_otu(scratch_path, 908125)

        new_otu_paths = set((scratch_path / "src").glob("*--*")) - pre_otu_paths

        assert len(new_otu_paths) == 1
        assert (new_otu_paths.pop() / "otu.json").exists()

    def test_failure(self, scratch_path: Path):
        """Test that adding an OTU with a pre-existing taxid fails."""
        pre_otu_paths = set((scratch_path / "src").glob("*--*"))

        self.run_add_otu(scratch_path, 345184)

        assert not set((scratch_path / "src").glob("*--*")) - pre_otu_paths


@pytest.mark.parametrize(
    "taxon_id,accession",
    [(908125, "NC_031754"), (1016856, "NC_015504")],
)
def test_init_and_add(taxon_id: int, accession: str, tmp_path: Path):
    path = tmp_path / "empty"

    init_reference(path)

    subprocess.run(
        [
            "virtool",
            "ref",
            "add",
            "otu",
            "--path",
            str(path),
            str(taxon_id),
        ],
        check=False,
    )

    build_path = tmp_path / "reference.json"

    build_json(False, build_path, path, "")

    subprocess.run(
        [
            "virtool",
            "ref",
            "build",
            "-o",
            str(build_path),
            "--path",
            str(path),
        ],
        check=True,
    )

    assert build_path.exists()
