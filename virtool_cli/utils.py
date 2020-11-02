import json
import os
from typing import Tuple

import aiofiles


def get_otu_paths(src_path: str) -> list:
    """
    Generates a list of paths to all OTUs in a src directory.

    :param src_path: Path to a src database directory
    :return: paths to all OTU in a src directory
    """
    alpha_paths = os.listdir(src_path)
    paths = []

    for alpha in alpha_paths:
        if alpha == "meta.json":
            continue

        otu_paths = [os.path.join(src_path, alpha, otu)
                     for otu in os.listdir(os.path.join(src_path, alpha))]

        for otu in otu_paths:
            paths.append(otu)

    return paths


async def get_otus(paths) -> dict:
    """
    Returns a mapping of every OTU path to their deserialized OTU dictionary.

    :param paths: List of paths to all OTU in a reference
    :return: A mapping of every OTU path to their taxid
    """
    taxids = dict()

    for otu_path in paths:
        async with aiofiles.open(os.path.join(otu_path, "otu.json"), "r") as f:
            otu = json.loads(await f.read())
            taxids[otu_path] = otu

    return taxids


async def get_isolates(path) -> dict:
    """
    Returns a mapping to every isolate and their folder name.

    :param path: A path to a OTU in a reference
    :return: A mapping of all of an OTU's isolates to their folder name (id)
    """
    isolates = dict()

    for folder in os.listdir(path):
        if folder != "otu.json" and "isolate.json" in set(os.listdir(os.path.join(path, folder))):
            async with aiofiles.open(os.path.join(path, folder, "isolate.json"), "r") as f:
                isolate = json.loads(await f.read())
                isolates[isolate["source_name"]] = folder

    return isolates


async def get_sequences(path) -> dict:
    """
    Returns a mapping of sequence accessions to their file name in a isolate directory.

    :param path: A path to an isolate in a reference
    :return: A mapping of all accessions in an isolate to their file name (id)
    """
    sequences = dict()

    for sequence_path in os.listdir(path):
        if sequence_path != "isolate.json":
            async with aiofiles.open(os.path.join(path, sequence_path), "r") as f:
                sequence = json.loads(await f.read())
                sequences[sequence["accession"]] = sequence_path

    return sequences


async def get_unique_ids(paths) -> Tuple[set, set]:
    """
    Returns sets containing unique random alphanumeric ids for both the isolates and the sequences

    :param paths: List of paths to all OTU in a reference
    :return: Sets containing unique ids for both isolates and sequences
    """
    isolate_ids = set()
    sequence_ids = set()

    for path in paths:
        for isolate_folder in os.listdir(path):
            if isolate_folder != "otu.json":
                isolate_ids.add(isolate_folder)
                for seq_file in os.listdir(os.path.join(path, isolate_folder)):
                    if seq_file != "isolate.json":
                        sequence_ids.add(seq_file.rstrip(".json"))

    return isolate_ids, sequence_ids
