import json
import os
from typing import Optional

from virtool_cli.utils import get_otu_paths, get_otus, create_otu_path


def run(src: str):
    """
    Fixes any folder-JSON name mismatches and incorrect taxid types

    :param src: Path to a given reference directory
    """
    paths = get_otu_paths(src)
    otus = get_otus(paths)
    otus_to_update = {}

    for path in paths:
        fix_folder_name(path, otus[path])
        otu = fix_taxid(otus[path])

        # if a new otu is returned then it should be updated
        if otu:
            otus_to_update[path] = otu

    write_otus(otus_to_update)


def fix_folder_name(path: str, otu: dict):
    """
    Fixes each OTU folder name by ensuring that it's the same as the "name" field in its otu.json file

    :param path: Path to a given reference directory
    :param otu: A deserialized otu.json
    """
    try:
        folder_name = path.rsplit("/", 1)[1]
        new_folder_name = create_otu_path(otu.get("name"))

        if folder_name != new_folder_name:
            new_path = path.replace(folder_name, new_folder_name)
            os.rename(path, new_path)
    except AttributeError:
        # OTU must not have name
        pass


def fix_taxid(otu: dict) -> Optional[dict]:
    """
    Ensures that each taxid inside every OTU's otu.json is of type int

    :param otu: A deserialized otu.json
    :return: The modified otu parameter if it needs to be updated, else None
    """
    try:
        taxid = otu["taxid"]
        if isinstance(taxid, str):
            return {
                **otu,
                "taxid": int(taxid)
            }
    except KeyError:
        # assure that taxid field is set to None
        return {
            **otu,
            "taxid": None
        }

    return None


def write_otus(otus: dict):
    """
    Write every otu.json file that's been modified with correct taxid type

    :param otus: A dictionary that maps file paths to each updated deserialized otu.json
    """
    for path in otus.keys():
        with open(os.path.join(path, "otu.json"), "w") as f:
            json.dump(otus[path], f, indent=4)
