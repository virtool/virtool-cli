import aiofiles
import json
import os


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


async def get_taxid_map(paths):
    """Returns a mapping of every OTU to their taxid."""
    taxids = {}
    for otu_path in paths:
        async with aiofiles.open(os.path.join(otu_path, "otu.json"), "r+") as f:
            otu = json.loads(await f.read())
            taxids[otu_path] = otu["taxid"]

    return taxids


async def get_isolates(path):
    """Returns a mapping to every isolate and their folder name."""
    isolates = {}
    for folder in os.listdir(path):
        if folder != "otu.json":
            if "isolate.json" in set(os.listdir(os.path.join(path, folder))):
                async with aiofiles.open(os.path.join(path, folder, "isolate.json"), "r") as f:
                    isolate = json.loads(await f.read())
                    isolates[isolate["source_name"]] = folder

    return isolates


async def get_sequences(path):
    """Returns a mapping of sequence accessions to their file name in a isolate directory."""
    sequences = {}
    for seq_file in os.listdir(path):
        if seq_file != "isolate.json":
            async with aiofiles.open(os.path.join(path, seq_file), "r") as f:
                sequence = json.loads(await f.read())
                sequences[sequence["accession"]] = seq_file

    return sequences


async def get_unique_ids(paths):
    """Returns a mapping of all isolate and accessions to their random alphanumeric id to ensure ids stay unique."""
    unique_ids = set()

    for path in paths:
        for isolate_folder in os.listdir(path):
            if isolate_folder != "otu.json":
                unique_ids.add(isolate_folder)
                for seq_file in os.listdir(os.path.join(path, isolate_folder)):
                    if seq_file != "isolate.json":
                        unique_ids.add(seq_file.rstrip(".json"))

    return unique_ids
