import shutil
import subprocess
from pathlib import Path

from syrupy import SnapshotAssertion

from virtool_cli.ref.legacy import Repo


def test_migrate(test_files_path: Path, snapshot: SnapshotAssertion, tmp_path: Path):
    """Test that the generated catalog creates all the same filenames as
    the control catalog
    """
    path = tmp_path / "v1"
    path.mkdir()

    migrate_path = path / "src"
    shutil.copytree(test_files_path / "src_v1", migrate_path)

    subprocess.call(["virtool", "ref", "migrate", "-src", str(migrate_path)])

    # Make sure the repo can be loaded.
    repo = Repo(path)

    # Make sure data from meta.json is available.
    assert repo.data_type == "genome"
    assert repo.organism == "virus"

    # Make sure the OTU directories are renamed.
    assert {
        otu_id: path.name for otu_id, path in repo.maps.otu_id_to_path.items()
    } == snapshot

    # Make sure exclusions.json files are created.
    assert all(
        (otu_path / "exclusions.json").exists()
        for otu_path in repo.maps.otu_id_to_path.values()
    )
