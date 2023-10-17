import pytest
import subprocess


@pytest.fixture()
def empty_repo(tmp_path):
    new_repo_path = tmp_path / "new_repo"

    subprocess.call(["virtool", "ref", "init", "-repo", str(new_repo_path)])

    return new_repo_path


def test_init_exist(empty_repo):
    assert empty_repo.exists()


def test_init_build(empty_repo):
    build_path = empty_repo / "reference.json"

    subprocess.call(
        [
            "virtool",
            "ref",
            "build",
            "-o",
            str(build_path),
            "-src",
            str(empty_repo / "src"),
        ]
    )

    assert build_path.exists()
