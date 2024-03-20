import asyncio
import json
import subprocess
from pathlib import Path

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from virtool_cli.add.accessions import add_accession
from virtool_cli.ref.build import build_json
from virtool_cli.ref.init import init_reference
from virtool_cli.utils.reference import get_sequence_paths


def get_all_sequence_paths(otu_path: Path) -> set[Path]:
    return {
        sequence_path
        for isolate_path in otu_path.iterdir()
        if isolate_path.is_dir()
        for sequence_path in get_sequence_paths(isolate_path)
    }


@pytest.mark.ncbi
class TestAddAccession:
    @staticmethod
    def get_all_filenames_in_dir(path: Path) -> set[str]:
        ignore_list = [".DS_Store", "Thumbs.db"]

        return {p.name for p in path.iterdir() if p.name not in ignore_list}

    def test_add_new_isolate(
        self,
        scratch_path: Path,
        snapshot: SnapshotAssertion,
    ):
        """Test that an accession is added when the parent isolate already exists."""
        otu_path = scratch_path / "src" / "cabbage_leaf_curl_jamaica_virus--d226290f"

        asyncio.run(add_accession("NC_038793", scratch_path))

        new_isolate_filenames = self.get_all_filenames_in_dir(otu_path) - {
            "496550f5",
            "d293d531",
            "exclusions.json",
            "otu.json",
        }

        assert len(new_isolate_filenames) == 1

        isolate_id = new_isolate_filenames.pop()

        with open(otu_path / isolate_id / "isolate.json") as f:
            assert json.load(f) == {
                "default": False,
                "id": isolate_id,
                "source_name": "CUc-3",
                "source_type": "isolate",
            }

        new_sequence_filenames = self.get_all_filenames_in_dir(
            otu_path / isolate_id
        ) - {
            "isolate.json",
        }

        assert len(new_sequence_filenames) == 1

        sequence_id = new_sequence_filenames.pop().replace(".json", "")

        with open(otu_path / isolate_id / f"{sequence_id}.json") as f:
            data = json.load(f)
            assert data == snapshot(name="sequence.json", exclude=props("_id"))
            assert data["_id"] == sequence_id
            assert data["accession"] == "NC_038793.1"

    def test_add_to_existing_isolate(
        self,
        scratch_path: Path,
        snapshot: SnapshotAssertion,
    ):
        """Test that an accession is added when the parent isolate already exists."""
        asyncio.run(add_accession("DQ178612", scratch_path))

        otu_path = scratch_path / "src" / "cabbage_leaf_curl_jamaica_virus--d226290f"

        assert not (
            self.get_all_filenames_in_dir(otu_path)
            - {
                "496550f5",
                "d293d531",
                "exclusions.json",
                "otu.json",
            }
        )

        new_sequence_filenames = self.get_all_filenames_in_dir(
            otu_path / "d293d531"
        ) - {
            "isolate.json",
            "ndaxyl7f.json",
            "pyf9x4wh.json",
        }

        assert len(new_sequence_filenames) == 1

        sequence_id = new_sequence_filenames.pop().replace(".json", "")

        with open(otu_path / "d293d531" / f"{sequence_id}.json") as f:
            data = json.load(f)
            assert data == snapshot(name="sequence.json", exclude=props("_id"))
            assert data["_id"] == sequence_id

    def test_accession_is_invalid(self, scratch_path: Path):
        """Test that adding an accession fails when it is invalid."""
        otu_path = (
            scratch_path / "src" / "pagoda_yellow_mosaic_associated_virus--dd21fd8f"
        )

        before = get_all_sequence_paths(otu_path)

        asyncio.run(add_accession("NC-024301", scratch_path))

        assert get_all_sequence_paths(otu_path) == before

    def test_accession_already_exists(
        self,
        scratch_path: Path,
    ):
        """Test that an accession cannot be added if it already exists."""
        otu_path = scratch_path / "src" / "abaca_bunchy_top_virus--c93ec9a9"

        before = get_all_sequence_paths(otu_path)

        asyncio.run(add_accession("NC_010319", scratch_path))

        assert get_all_sequence_paths(otu_path) == before


@pytest.mark.ncbi
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


@pytest.mark.ncbi
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


@pytest.mark.ncbi
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
