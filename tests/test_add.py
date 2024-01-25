import subprocess
from pathlib import Path

import pytest

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
    def run_command(accession: str, src_path: Path):
        subprocess.run(
            [
                "virtool",
                "ref",
                "add",
                "accession",
                "-acc",
                accession,
                "-src",
                str(src_path),
            ],
            check=False,
        )

    def run_add_accession(self, accession: str, otu_dirname: str, src_path: Path):
        """Add into an existing isolate directory"""
        otu_path = src_path / otu_dirname

        pre_sequence_paths = get_all_sequence_paths(otu_path)

        self.run_command(accession, src_path)

        post_sequence_paths = get_all_sequence_paths(otu_path)

        return pre_sequence_paths, post_sequence_paths

    def run_add_isolate(self, accession: str, otu_subpath: str, src_path: Path):
        """Add into a new isolate directory"""
        otu_path = src_path / otu_subpath

        pre_isolate_paths = set(get_isolate_paths(otu_path))

        self.run_command(accession, src_path)

        post_isolate_paths = set(get_isolate_paths(otu_path))

        return pre_isolate_paths, post_isolate_paths

    @pytest.mark.parametrize(
        "accession, otu_subpath",
        [
            ("DQ178612", "cabbage_leaf_curl_jamaica_virus--d226290f"),
            ("NC_038793", "cabbage_leaf_curl_jamaica_virus--d226290f"),
        ],
    )
    def test_add_accession_success(
        self,
        accession: str,
        otu_subpath,
        src_scratch_path: Path,
    ):
        """Check that virtool ref add accession does the job when the isolate exists."""
        pre_sequence_paths, post_sequence_paths = self.run_add_accession(
            accession,
            otu_dirname=otu_subpath,
            src_path=src_scratch_path,
        )

        assert post_sequence_paths.difference(pre_sequence_paths)

    @pytest.mark.parametrize(
        "accession, isolate_subpath",
        [
            ("NC_010319", "abaca_bunchy_top_virus--c93ec9a9/4e8amg20"),
            ("NC-024301", "pagoda_yellow_mosaic_associated_virus--dd21fd8f"),
        ],
    )
    def test_add_accession_fail(
        self,
        accession,
        isolate_subpath,
        src_scratch_path: Path,
    ):
        """Check that virtool ref add accession does not add sequences that already exist"""
        pre_sequence_paths, post_sequence_paths = self.run_add_accession(
            accession,
            otu_dirname=isolate_subpath,
            src_path=src_scratch_path,
        )

        new_sequences = post_sequence_paths.difference(pre_sequence_paths)

        assert not new_sequences

    @pytest.mark.parametrize(
        "accession, otu_dirname",
        [("KT390494", "nanovirus_like_particle--ae0f2a35")],
    )
    def test_add_isolate_success(self, accession, otu_dirname, src_scratch_path: Path):
        """Check that virtool ref add accession does the job when the isolate does not exist"""
        otu_path = src_scratch_path / otu_dirname

        pre_isolate_paths = set(get_isolate_paths(otu_path))

        self.run_command(accession, src_scratch_path)

        post_isolate_paths = set(get_isolate_paths(otu_path))

        new_isolates = post_isolate_paths.difference(pre_isolate_paths)

        assert new_isolates
        assert (new_isolates.pop() / "isolate.json").exists()


class TestAddAccessions:
    @staticmethod
    def run_command(accessions: str, otu_path: Path):
        subprocess.run(
            [
                "virtool",
                "ref",
                "add",
                "accessions",
                "--accessions",
                accessions,
                "--otu_path",
                str(otu_path),
            ],
            check=False,
        )

    def run_add_accessions(self, accessions: str, otu_dirname: str, src_path: Path):
        """Add into an existing isolate directory"""
        otu_path = src_path / otu_dirname

        pre_sequence_paths = get_all_sequence_paths(otu_path)

        self.run_command(accessions, otu_path)

        post_sequence_paths = get_all_sequence_paths(otu_path)

        return pre_sequence_paths, post_sequence_paths

    @pytest.mark.parametrize(
        "accessions, otu_subpath",
        [
            ("DQ178612, NC_038793", "cabbage_leaf_curl_jamaica_virus--d226290f"),
            ("KT390494, KT390496, KT390501", "nanovirus_like_particle--ae0f2a35"),
        ],
    )
    def test_add_accessions_success(
        self,
        accessions,
        otu_subpath,
        src_scratch_path: Path,
    ):
        """Check that virtool ref add accessions does the job"""
        pre_sequence_paths, post_sequence_paths = self.run_add_accessions(
            accessions,
            otu_dirname=otu_subpath,
            src_path=src_scratch_path,
        )

        new_sequences = post_sequence_paths.difference(pre_sequence_paths)

        assert new_sequences

    @pytest.mark.parametrize(
        "accessions, otu_subpath",
        [
            ("DQ178611, DQ178610", "cabbage_leaf_curl_jamaica_virus--d226290f"),
            ("KTX390494, KTX390496, KTX390501", "nanovirus_like_particle--ae0f2a35"),
        ],
    )
    def test_add_accessions_fail(self, accessions, otu_subpath, src_scratch_path: Path):
        """Check that virtool ref add accessions does the job"""
        pre_sequence_paths, post_sequence_paths = self.run_add_accessions(
            accessions,
            otu_dirname=otu_subpath,
            src_path=src_scratch_path,
        )

        assert post_sequence_paths.union(pre_sequence_paths) == pre_sequence_paths


class TestAddOTU:
    @staticmethod
    def run_command(taxon_id: int, src_path: Path):
        subprocess.run(
            [
                "virtool",
                "ref",
                "add",
                "otu",
                "-taxid",
                str(taxon_id),
                "-src",
                str(src_path),
            ],
            check=False,
        )

    def run_add_otu(self, taxon_id: int, src_path: Path):
        """Attempt to add a new OTU"""
        pre_otu_paths = set(src_path.glob("*--*"))

        self.run_command(taxon_id, src_path)

        post_otu_paths = set(src_path.glob("*--*"))

        return pre_otu_paths, post_otu_paths

    def test_add_otu_success(self, src_scratch_path: Path):
        """ """
        taxon_id = 908125

        pre_otu_paths = set(src_scratch_path.glob("*--*"))

        self.run_add_otu(taxon_id=taxon_id, src_path=src_scratch_path)

        post_otu_paths = set(src_scratch_path.glob("*--*"))

        new_otus = post_otu_paths.difference(pre_otu_paths)

        assert new_otus

        otu_path = new_otus.pop()

        assert (otu_path / "otu.json").exists()

    @pytest.mark.parametrize("taxon_id", [345184])
    def test_add_otu_conflict(self, taxon_id, src_scratch_path: Path):
        """Don't add taxon IDs that already exist"""
        pre_otu_paths = set(src_scratch_path.glob("*--*"))

        self.run_add_otu(taxon_id=taxon_id, src_path=src_scratch_path)

        post_otu_paths = set(src_scratch_path.glob("*--*"))

        new_otus = post_otu_paths.difference(pre_otu_paths)

        assert not new_otus


@pytest.mark.parametrize(
    "taxon_id, accession",
    [(908125, "NC_031754"), (1016856, "NC_015504")],
)
def test_init_and_add(taxon_id: int, accession: str, tmp_path: Path):
    empty_repo_path = tmp_path / "empty"

    init_reference(empty_repo_path)

    src_path = empty_repo_path / "src"

    TestAddOTU.run_command(taxon_id, src_path)

    TestAddAccessions.run_command(accession, src_path)

    build_path = tmp_path / "reference.json"

    subprocess.run(
        [
            "virtool",
            "ref",
            "build",
            "-o",
            str(build_path),
            "-src",
            str(src_path),
        ],
        check=False,
    )

    assert build_path.exists()
