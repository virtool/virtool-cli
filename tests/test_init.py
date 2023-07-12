import pytest
import os
from pathlib import Path
from dotenv import load_dotenv
from requests import JSONDecodeError

from paths import TEST_FILES_PATH

TEST_PATH = TEST_FILES_PATH / "new_repo"
env = load_dotenv(TEST_FILES_PATH / 'testing.env')

cache_path = TEST_FILES_PATH / '.cache'
cache_path.mkdir(exist_ok=True)

from virtool_cli import init

def test_bad_extant_url():
    with pytest.raises(KeyError):
        caller_tarfile = init.fetch_callers(
            api_url='https://api.github.com',
            cache_path=cache_path)

if __name__ == '__main__':
    print(os.getenv("CALLER_API_URL"))