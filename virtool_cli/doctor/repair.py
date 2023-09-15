import json
from typing import Optional
from pathlib import Path
# import structlog
import logging
import re

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import get_otu_paths, generate_otu_dirname


def run(src_path: Path, debugging: bool = False):
    """
    Fixes any folder-JSON name mismatches and incorrect taxid types

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    repair_reference(src_path)
    
def repair_reference(src_path):
    """
    Runs legacy reference repair functions

    :param src_path: Path to a given reference directory
    """
    otu_paths = get_otu_paths(src_path)
    # Rename OTU folders
    otus = map_otus(otu_paths)
    otus_to_update = {}
    for otu_path in otu_paths:
        results = []

        new_path = fix_folder_name(otu_path, otus.get(otu_path))
        # if a folder name has been changed then a new path will return
        if new_path:
            results.append("Fixed folder name")
            # making sure to update the dictionary that maps paths to OTUs
            otus[new_path] = otus.get(otu_path)
            otu_path = new_path

        otu = fix_taxid(otus[otu_path])
        # if a new otu is returned then it should be updated
        if otu:
            otus_to_update[otu_path] = otu
            results.append("Fixed taxid field")
            
        log_results(results, otus[otu_path]["name"])
    write_otus(otus_to_update)


def fix_folder_name(path: Path, otu: dict) -> Optional[str]:
    """
    Fixes each OTU folder name by ensuring that it's the same as the "name" field in its otu.json file

    :param path: Path to a given reference directory
    :param otu: A deserialized otu.json
    """
    new_folder_name = generate_otu_dirname(otu['name'], otu['_id'])

    if path.name != new_folder_name:
        new_path = path.with_name(new_folder_name)
        path.rename(new_path)
        return new_path


def fix_taxid(otu: dict) -> Optional[dict]:
    """
    Ensures that each taxid inside every OTU's otu.json is of type int

    :param otu: A deserialized otu.json
    :return: The modified otu parameter if it needs to be updated, else None
    """
    try:
        taxid = otu.get("taxid", None)
        if isinstance(taxid, str):
            return {**otu, "taxid": int(taxid)}
    except KeyError:
        # assure that taxid field is set to None
        return {**otu, "taxid": None}

    return None


def log_results(results: list, name: str):
    """
    Log repair results to console

    :param name: Name of OTU to log results for
    :param results: List of successful operations performed for a given OTU
    :param console: Rich console object used to log results
    :return:
    """
    if not results:
        return

    base_logger.info(f"Changed {name}", name=name, result=results)


def write_otus(otus: dict):
    """
    Write every otu.json file that's been modified with correct taxid type

    :param otus: A dictionary that maps file paths to each updated deserialized otu.json
    """
    for path in otus.keys():
        with open(path / "otu.json", "w") as f:
            json.dump(otus.get(path), f, indent=4)

def map_otus(paths: list) -> dict:
    """
    Returns a mapping of every OTU path to their deserialized OTU dictionary.

    :param paths: List of paths to all OTU in a reference
    :return: A mapping of every OTU path to its OTU dictionary
    """
    path_taxid_map = {}

    for path in paths:
        path_taxid_map[path] = json.loads((path / "otu.json").read_text())

    return path_taxid_map