import json
from typing import Optional
import pathlib

from rich.console import Console

from virtool_cli.utils import get_otu_paths, get_otus, create_otu_path


def run(src_path: pathlib.Path):
    """
    Fixes any folder-JSON name mismatches and incorrect taxid types

    :param src_path: Path to a given reference directory
    """
    paths = get_otu_paths(src_path)
    otus = get_otus(paths)
    otus_to_update = {}
    console = Console()

    for path in paths:
        results = []

        new_path = fix_folder_name(path, otus[path])
        # if a folder name has been changed then a new path will return
        if new_path:
            results.append("    - Fixed folder name")
            # making sure to update the dictionary that maps paths to OTUs
            otus[new_path] = otus[path]
            path = new_path

        otu = fix_taxid(otus[path])
        # if a new otu is returned then it should be updated
        if otu:
            otus_to_update[path] = otu
            results.append("    - Fixed taxid field")

        log_results(results, otus[path]["name"], console)

    write_otus(otus_to_update)


def fix_folder_name(path: pathlib.Path, otu: dict) -> Optional[str]:
    """
    Fixes each OTU folder name by ensuring that it's the same as the "name" field in its otu.json file

    :param path: Path to a given reference directory
    :param otu: A deserialized otu.json
    """
    new_folder_name = create_otu_path(otu.get("name"))

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
        taxid = otu["taxid"]
        if isinstance(taxid, str):
            return {**otu, "taxid": int(taxid)}
    except KeyError:
        # assure that taxid field is set to None
        return {**otu, "taxid": None}

    return None


def log_results(results: list, name: str, console: Console):
    """
    Log repair results to console

    :param name: Name of OTU to log results for
    :param results: List of successful operations performed for a given OTU
    :param console: Rich console object used to log results
    :return:
    """
    if not results:
        return

    console.print(f"[green]✔ {name}")

    for result in results:
        console.print(f"[green]{result}")
    console.print()


def write_otus(otus: dict):
    """
    Write every otu.json file that's been modified with correct taxid type

    :param otus: A dictionary that maps file paths to each updated deserialized otu.json
    """
    for path in otus.keys():
        with open(path / "otu.json", "w") as f:
            json.dump(otus[path], f, indent=4)
