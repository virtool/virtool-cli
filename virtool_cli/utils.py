import aiofiles
import json
import os


def get_paths(src_path: str) -> list:
    """
    Generates a list of paths to all OTU in a src directory.

    Parameters:
        src_path (str): Path to a src database directory

    Returns:
        A list containing paths to all OTU in a src directory

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


async def get_taxids(paths):
    """Returns a mapping of every OTU to their taxid"""
    taxids = {}
    for otu_path in paths:
        async with aiofiles.open(os.path.join(otu_path, "otu.json"), "r+") as f:
            otu = json.loads(await f.read())
            taxids[otu_path] = otu["taxid"]

    return taxids


async def get_isolates(path):
    """Returns a mapping to every isolate and their folder name"""
    isolates = {}
    for folder in os.listdir(path):
        if folder != "otu.json":
            if "isolate.json" in set(os.listdir(os.path.join(path, folder))):
                async with aiofiles.open(os.path.join(path, folder, "isolate.json"), "r") as f:
                    isolate = json.loads(await f.read())
                    isolates[isolate["source_name"]] = folder

    return isolates
