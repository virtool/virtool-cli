import pytest
import os
import json
import subprocess

from paths import TEST_FILES_PATH

TEST_SRC_PATH = TEST_FILES_PATH / "src_a"
TEST_ACCLOG_PATH = TEST_FILES_PATH / "catalog"


@pytest.fixture(scope="session", autouse=True)
def command():
    command = [
        "virtool", "ref", "update", 
        "-src", TEST_SRC_PATH, 
        "-cat", TEST_ACCLOG_PATH
    ]
    subprocess.call(command)
