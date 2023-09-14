import pytest
import shutil
import json
import subprocess
from pathlib import Path

from paths import TEST_FILES_PATH
BASE_PATH = TEST_FILES_PATH / "src_test_partial"
TEST_ACCLOG_PATH = TEST_FILES_PATH / "catalog"


def convert_to_dict(otu_list: list) -> dict:
    """
    Converts a list of OTUs to a dict keyed by Virtool ID
    
    :param otu_list: A list of deserialized OTU data
    :return: The contents of otu_list keyed by OTU id 
    """
    otu_dict = {}
    for otu in otu_list:
        otu_dict[otu['_id']] = otu
    return otu_dict

def get_otu_accessions(otu_dict: dict) -> set:
    """
    Gets all accessions from an OTU directory and returns a list

    :param otu_dict: Deserialized OTU data
    :return: The accessions included under the OTU directory in a set
    """
    accessions = set()
    
    for isolate in otu_dict['isolates']:
        for sequence in isolate['sequences']:
            accessions.add(sequence['accession'])

    return accessions

@pytest.mark.parametrize("base_path", [BASE_PATH])
def test_update(base_path, tmp_path):
    """
    Test that updates actually pull something
    """
    fetch_path = tmp_path / 'src'
    pre_update_ref_path = tmp_path / 'reference_pre.json'
    post_update_ref_path = tmp_path / 'reference_post.json'

    shutil.copytree(base_path, fetch_path)

    subprocess.call([
        'virtool', 'ref', 'build', 
        '-src', str(fetch_path), 
        '-o', str(pre_update_ref_path)
    ])

    subprocess.call([
        "virtool", "ref", "update", "reference",
        "-src", str(fetch_path),
        "-cat", str(TEST_ACCLOG_PATH),
    ])

    subprocess.call([
        'virtool', 'ref', 'build', 
        '-src', str(fetch_path), 
        '-o', str(post_update_ref_path)
    ])
    
    
    reference_pre = json.loads(pre_update_ref_path.read_text())
    pre_otu_dict = convert_to_dict(reference_pre['otus'])

    reference_post = json.loads(post_update_ref_path.read_text())
    post_otu_dict = convert_to_dict(reference_post['otus'])

    difference_counter = 0

    for otu_id in post_otu_dict:
        pre_accessions = get_otu_accessions(pre_otu_dict[otu_id])
        post_accessions = get_otu_accessions(post_otu_dict[otu_id])

        print(pre_accessions)
        print(post_accessions)

        if pre_accessions != post_accessions:
             difference_counter += 1

    # Any new data counts
    assert difference_counter > 0