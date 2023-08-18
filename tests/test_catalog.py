import pytest
from pathlib import Path
import shutil
import json
import subprocess
import logging

from paths import TEST_FILES_PATH
TEST_FULLSRC_PATH = TEST_FILES_PATH / "src_test"
TEST_PARTSRC_PATH = TEST_FILES_PATH / "src_test_partial"
TEST_CONTROL_LOG_PATH = TEST_FILES_PATH / "catalog"

CONTROL_SET = { child.name for child in Path(TEST_CONTROL_LOG_PATH).iterdir() }

LOGGER = logging.getLogger(__name__)

@pytest.fixture()
def catalog_path(tmp_path):
    return tmp_path / "catalog"

@pytest.fixture()
def command(src_path, catalog_path):
    return [
        "virtool", "acc", "init", 
        "-src", str(src_path),
        "-cat", str(catalog_path),
    ]

# @pytest.mark.parametrize("src_path", [TEST_FULLSRC_PATH])
# def test_init_full(command, src_path, catalog_path):
#     """
#     Test that the generated catalog creates all the same filenames as
#     the control catalog
#     """

#     subprocess.call(command)

#     temp_catalog = catalog_path

#     temp_set = { child.name for child in temp_catalog.iterdir() }

#     LOGGER.info(f'generated: /n{temp_set}')

#     LOGGER.info(f'control: /n {CONTROL_SET}')

#     assert temp_set == CONTROL_SET

#     for otu_path in src_path.iterdir():
#         if not otu_path.is_dir(): continue
        
#         [ _, otu_id ] = otu_path.name.split('--')
#         assert len([temp_catalog.glob(f'*--{otu_id}.json')]) > 0        

@pytest.mark.parametrize("src_path", [TEST_PARTSRC_PATH])
def test_init_partial(command, src_path, catalog_path):
    """
    Test that the generated catalog generates a subset of 
    the control catalog
    """
    subprocess.call(command)

    temp_catalog = catalog_path

    temp_set = { child.name for child in temp_catalog.iterdir() }

    LOGGER.info(f'generated: /n{temp_set}')

    LOGGER.info(f'control: /n {CONTROL_SET}')

    assert temp_set.issubset(CONTROL_SET)

    for otu_path in src_path.iterdir():
        if not otu_path.is_dir(): continue
        
        [ _, otu_id ] = otu_path.name.split('--')
        assert len([temp_catalog.glob(f'*--{otu_id}.json')]) > 0

@pytest.mark.parametrize("src_path", [TEST_FULLSRC_PATH])
def test_update(src_path, catalog_path):
    """
    Check that update can replace missing listings
    """
    shutil.copytree(TEST_CONTROL_LOG_PATH, catalog_path)
    first_listing = list(catalog_path.glob('*.json'))[0]
    try:
        first_listing.unlink()
    except Exception as e:
        print(e)

    subprocess.call([
        "virtool", "acc", "update", 
        "-src", str(src_path),
        "-cat", str(catalog_path),
    ])

    updated_listings = [ path.name for path in catalog_path.glob('*.json') ]

    for path in TEST_CONTROL_LOG_PATH.glob('*.json'):
        assert path.name in updated_listings

@pytest.mark.parametrize("src_path", [TEST_PARTSRC_PATH])
def test_check_contents(command, catalog_path):
    """
    """
    subprocess.call(command)

    temp_catalog = catalog_path

    for path in temp_catalog.glob('*.json'):
        with open(TEST_CONTROL_LOG_PATH / path.name, "r") as f:
            control_listing_data = json.load(f)

        gen_listing_data = json.loads(path.read_text())

        gen_listing_data.pop('accessions')
        control_listing_data.pop('accessions')

        assert gen_listing_data == control_listing_data