import pytest
import filecmp


from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_curate import *
from virtool_cli.vfam_polyprotein import *
from virtool_cli.vfam_markov import *
from pathlib import Path


LARGE_INPUT_PATH = "large_input"
ABC_FILE_PATH = "vfam_og_intermediates/abc_file.abc"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])