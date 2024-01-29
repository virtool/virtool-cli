from pathlib import Path

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from virtool_cli.repo.cls import Repo


def test_init(scratch_path: Path, snapshot: SnapshotAssertion):
    """Test that the repo is initialized correctly."""
    repo = Repo(scratch_path)

    assert repo.data_type == "genome"
    assert repo.organism == "virus"
    assert repo.path == scratch_path

    assert {
        key: path.name for key, path in repo.maps.otu_id_to_path.items()
    } == snapshot(name="otu_id_to_path")

    assert repo.maps.sequence_id_to_otu_id == snapshot(name="sequence_id_to_otu_id")

    assert {
        key: path.name for key, path in repo.maps.sequence_id_to_path.items()
    } == snapshot(name="sequence_id_to_path")

    assert repo.maps.taxid_to_otu_id == snapshot(name="taxid_to_otu_id")


class TestGetOTUByID:
    def test_ok(self, scratch_path: Path, snapshot: SnapshotAssertion):
        """Test that an OTU is retrieved correctly."""
        otu = Repo(scratch_path).get_otu_by_id("3962f6ec")

        assert otu == snapshot(name="otu", exclude=props("path"))

        assert otu.path.name == "jacquemontia_mosaic_yucatan_virus--3962f6ec"
        assert otu.path.parent == scratch_path / "src"

    def test_not_found(self, scratch_path: Path):
        """Test that a ValueError is raised when the OTU is not found."""
        with pytest.raises(ValueError) as e:
            Repo(scratch_path).get_otu_by_id("foobar")

        assert str(e.value) == "No OTU with ID foobar"


class TestGetOTUByTaxid:
    def test_ok(self, scratch_path: Path, snapshot: SnapshotAssertion):
        """Test that an OTU is retrieved correctly."""
        otu = Repo(scratch_path).get_otu_by_taxid(1278205)

        assert otu == snapshot(name="otu", exclude=props("path"))

        assert otu.path.name == "dahlia_latent_viroid--ab6bd88a"
        assert otu.path.parent == scratch_path / "src"

    def test_not_found(self, scratch_path: Path):
        """Test that a ValueError is raised when the OTU is not found."""
        with pytest.raises(ValueError) as e:
            Repo(scratch_path).get_otu_by_taxid(512)

        assert str(e.value) == "No OTU with taxid 512"


def test_create_otu(scratch_path: Path, snapshot: SnapshotAssertion):
    """Test that an OTU is created in the correct location."""
    repo = Repo(scratch_path)

    otu = repo.create_otu("Tobacco mosaic virus", 512, "TMV")

    assert otu == snapshot(name="otu", exclude=props("id", "path"))

    assert otu.path.name == f"tobacco_mosaic_virus--{otu.id}"
    assert otu.path.parent.name == "src"
    assert otu.path.parent.parent == repo.path == scratch_path

    assert repo.maps.otu_id_to_path[otu.id] == otu.path
    assert repo.maps.taxid_to_otu_id[512] == otu.id
