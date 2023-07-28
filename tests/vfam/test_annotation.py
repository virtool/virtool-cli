import pytest
import json

from fixtures import *
from paths_vfam import VFAM_INTERMEDIATES_PATH

DUPES_JSON = VFAM_INTERMEDIATES_PATH / "Dupes" / "master.json"
GENERIC_JSON = VFAM_INTERMEDIATES_PATH / "Generic" / "master.json"
LARGE_JSON = VFAM_INTERMEDIATES_PATH / "Large" / "master.json"


@pytest.mark.parametrize("input_path", [DUPES_JSON, GENERIC_JSON, LARGE_JSON])
def test_get_taxonomy(input_path):

    with input_path.open("r") as handle:
        json_decode = json.load(handle)
        for annotation in json_decode:
            family_count = 0
            if "families" in annotation:
                for i in annotation["families"]:
                    family_count += annotation["families"][i]
                assert family_count == len(annotation["entries"])
                assert annotation["count"] == family_count

            genus_count = 0
            if "genera" in annotation:
                for i in annotation["genera"]:
                    genus_count += annotation["genera"][i]
                assert genus_count == family_count
