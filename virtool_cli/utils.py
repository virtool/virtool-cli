import json
import os
from typing import Tuple
import pathlib

import aiofiles
from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

NCBI_REQUEST_INTERVAL = 0.3 if Entrez.email and Entrez.api_key else 0.6


def get_otu_paths(src_path: pathlib.Path) -> list:
    """
    Generates a list of paths to all OTUs in a src directory.

    :param src_path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    paths = []

    for alpha in src_path.iterdir():
        if alpha.name == "meta.json":
            continue

        otu_paths = [otu for otu in alpha.iterdir()]
        for otu in otu_paths:
            paths.append(otu)

    return paths


def get_otus(paths: list) -> dict:
    """
    Returns a mapping of every OTU path to their deserialized OTU dictionary.

    :param paths: List of paths to all OTU in a reference
    :return: A mapping of every OTU path to its OTU dictionary
    """
    path_taxid_map = dict()

    for path in paths:
        with open(path / "otu.json", "r") as f:
            otu = json.load(f)
            path_taxid_map[path] = otu

    return path_taxid_map


def create_otu_path(
    otu_name: str, reference_path: pathlib.Path = None, first_letter: str = None
) -> pathlib.Path:
    """
    Generates a new path in a reference for an OTU directory if the full path is provided, else it just returns
    the formatted name of an OTU directory

    :param otu_name: Lowercase name of an OTU to be appended on the end of the path
    :param reference_path: Path to a reference directory
    :param first_letter: First letter of an OTU name
    :return: Path in a reference to generate an OTU directory
    """
    if reference_path and first_letter:
        return (
            reference_path
            / first_letter
            / otu_name.replace(" ", "_").replace("/", "_").lower()
        )

    return otu_name.replace(" ", "_").replace("/", "_").lower()


async def get_isolates(path: pathlib.Path) -> dict:
    """
    Returns a mapping to every isolate and their folder name.

    :param path: A path to a OTU directory in a reference
    :return: A mapping of all of an OTU's isolates to their folder name (id)
    """
    isolates = dict()

    for folder in path.iterdir():
        # ignore the otu.json file in the OTU folder and parse through isolate folders
        if folder.name != "otu.json" and folder / "isolate.json" in folder.iterdir():
            async with aiofiles.open(folder / "isolate.json", "r") as f:
                isolate = json.loads(await f.read())
                isolates[isolate["source_name"]] = folder.name

    return isolates


async def get_sequences(path: pathlib.Path) -> dict:
    """
    Returns a mapping of sequence accessions to their file name in a isolate directory.

    :param path: A path to an isolate directory in a reference
    :return: A mapping of all accessions in an isolate to their file name (id)
    """
    sequences = dict()

    for sequence_id in path.iterdir():
        if sequence_id.name != "isolate.json":
            async with aiofiles.open(sequence_id, "r") as f:
                sequence = json.loads(await f.read())
                sequences[sequence["accession"]] = sequence_id.name

    return sequences


async def get_unique_ids(paths: list) -> Tuple[set, set]:
    """
    Returns sets containing unique random alphanumeric ids for both the isolates and the sequences

    :param paths: List of paths to all OTU in a reference
    :return: Sets containing unique ids for both isolates and sequences
    """
    isolate_ids = set()
    sequence_ids = set()

    for path in paths:
        for isolate_id in path.iterdir():
            if isolate_id.name != "otu.json":
                isolate_ids.add(isolate_id)
                for seq_id in isolate_id.iterdir():
                    if seq_id.name != "isolate.json":
                        sequence_ids.add(seq_id.name.rstrip(".json"))

    return isolate_ids, sequence_ids
