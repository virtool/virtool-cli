import pytest
from pathlib import Path
import json
import shutil
import subprocess

from paths import TEST_FILES_PATH
SRC_V1_PATH = TEST_FILES_PATH / "src_v1"
REF_V1_PATH = TEST_FILES_PATH / "reference_v1.json"


def convert_to_dict(otu_list: list):
    otu_dict = {}
    for otu in otu_list:
        otu_dict[otu['_id']] = otu
    return otu_dict

@pytest.mark.parametrize("base_path", [SRC_V1_PATH])
def test_migrate(base_path, tmp_path):
    """
    Test that the generated catalog creates all the same filenames as
    the control catalog
    """
    migrate_path = tmp_path / 'migration'

    shutil.copytree(base_path, migrate_path)

    assert len([migrate_path.glob('[a-z]')]) > 0

    reference_v1 = json.loads(REF_V1_PATH.read_text())

    v1_otu_set = convert_to_dict(reference_v1['otus'])

    ref_v2_path = tmp_path / 'reference_v2.json'

    subprocess.call(['virtool', 'ref', 'migrate', '-src', str(migrate_path) ])

    subprocess.call([
        'virtool', 'ref', 'build', 
        '-src', str(migrate_path), 
        '-o', str(ref_v2_path)
    ])

    reference_v2 = json.loads(ref_v2_path.read_text())

    v2_otu_set = convert_to_dict(reference_v2['otus'])

    for otu_id in v2_otu_set:
        for key in v2_otu_set[otu_id]:
            assert v2_otu_set[otu_id][key] == v1_otu_set[otu_id][key]