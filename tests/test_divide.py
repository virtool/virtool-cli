import filecmp
import subprocess
from pathlib import Path

import pytest


@pytest.mark.parametrize("src", ["reference.json", "reference_with_indent.json"])
def test_divide(src, src_test_path: Path, test_files_path: Path, tmpdir):
    """
    Tests the divide operation to see if it produces a valid directory

    """
    output_path = tmpdir / "divide_output"

    subprocess.call(
        [
            "virtool",
            "ref",
            "divide",
            "-o",
            str(output_path),
            "-f",
            str(test_files_path / src),
        ]
    )

    cmp = filecmp.dircmp(src_test_path, output_path)

    assert len(cmp.diff_files) == 0
