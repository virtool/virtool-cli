import asyncio
import json
import os
import time
from concurrent.futures.thread import ThreadPoolExecutor
from random import choice
from string import ascii_letters, ascii_lowercase, digits
from typing import Union, Iterable, Tuple
from urllib.error import HTTPError

import aiofiles
import aiojobs
from Bio import Entrez, SeqIO
from rich.console import Console

from virtool_cli.utils import get_otu_paths, get_taxid_map, get_isolates, get_sequences, get_unique_ids

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

REQUEST_INTERVAL = 0.4 if Entrez.email and Entrez.api_key else 0.6


async def isolate(src):
    """
    Runs routines to find new isolates for OTU in a reference directory and writes newfound accessions to local cache

    :param src: Path to a given reference directory
    """
    scheduler = await aiojobs.create_scheduler(limit=10)
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor())

    q = asyncio.Queue()

    existing_accessions = get_cache()

    # get paths for all OTU in the directory
    paths = get_otu_paths(src)

    # get mapping of all OTU paths to their taxid
    taxid_otu_path_map = await get_taxid_map(paths)

    for path in paths:
        taxid = str(taxid_otu_path_map[path])

        if taxid is None:
            continue

        accessions = existing_accessions.get(taxid)

        await scheduler.spawn(fetch_otu_isolates(taxid, path, accessions, paths, q))

        await asyncio.sleep(REQUEST_INTERVAL)

    results_to_cache = list()
    while True:
        try:
            results_to_cache.append(q.get_nowait())
        except asyncio.QueueEmpty:
            if not scheduler.active_count:
                break

        await asyncio.sleep(0.01)

    await scheduler.close()

    updated_cache = update_cache(existing_accessions, results_to_cache)

    write_cache(updated_cache)


async def fetch_otu_isolates(taxid, path, accessions, paths, queue):
    """
    Routine to fetch new isolates for a given OTU, creating a directory for each

    :param taxid: The taxon id for a given OTU
    :param path: Path to a given OTU in a reference directory
    :param accessions: List of accessions for the given taxon id from local cache
    :param paths: List of all paths to every OTU in a reference
    :param queue: Asyncio queue to put newfound accessions onto for caching purposes
    """

    isolates = await get_isolates(path)

    records, new_accessions = await asyncio.get_event_loop().run_in_executor(None, get_records, accessions, taxid)

    isolate_ids, sequence_ids = await get_unique_ids(paths)

    if records is not None:
        for accession in [accession for accession in records.values() if accession.seq]:
            accession_data = await get_qualifiers(accession.features)

            # try to find isolate type and name automatically
            isolate_type = await find_isolate(accession_data)

            if isolate_type is None:
                continue

            # see if isolate directory already exists and create it if needed

            isolate_name = accession_data.get(isolate_type)[0]

            if isolate_name not in isolates:
                new_id = await store_isolate(path, isolate_name, isolate_type, isolate_ids)
                isolate_ids.add(new_id)
                isolates[isolate_name] = new_id

            # create new sequence file
            isolate_path = os.path.join(path, isolates.get(isolate_name))
            new_id = await store_sequence(isolate_path, accession, accession_data, sequence_ids)
            sequence_ids.add(new_id)

            await queue.put((taxid, new_accessions))


def update_cache(existing_cache, results) -> dict:
    """
    Updates the existing accessions cache with fetch results

    :param existing_cache: The existing cache from the .cli folder
    :param results: List of tuples pulled out of the asyncio Queue containing each taxid and their accessions
    """
    new_accessions = dict()
    for taxid, accessions in results:
        new_accessions[taxid] = accessions

    existing_cache.update(new_accessions)

    return existing_cache


def get_records(accessions, taxid) -> Union[Tuple[dict, list], Tuple[None, None]]:
    """
    Search the NCBI database for accession ids if not found in cache, then fetches all accession records

    :param accessions: A list of accession ids if found in cache, else None
    :param taxid: Taxon id for a given OTU
    :return: A dictionary with all the accession records and a list containing all accession ids, or None if records aren't found
    """
    console = Console()

    if accessions is None:
        accessions = []

        while True:
            try:
                record = Entrez.read(Entrez.elink(dbfrom="taxonomy", db="nucleotide",
                                                  id=taxid, idtype="acc"))
                break
            except HTTPError as e:
                if e.code == 429:
                    time.sleep(2)

        for linksetdb in record[0]["LinkSetDb"][0]["Link"]:
            accessions.append(linksetdb["Id"])

    while True:
        try:
            records = fetch_records(accessions)
            console.print(f"Found isolate data for {taxid}", style="green")
            break
        except HTTPError as e:
            if e.code == 429:
                time.sleep(2)
            elif e.code == 500:
                console.print(f"Could not find isolate data for {taxid}", style="red")
                return None, None

    return records, accessions


def fetch_records(accessions) -> dict:
    """
    Fetch all accession ids given by a list, then converts the records to a SeqIO dictionary object

    :param accessions: List of accession ids
    :return: A SeqIO dictionary object
    """
    handle = Entrez.efetch(db="nucleotide", id=accessions,
                           rettype="gb", retmode="text")

    return SeqIO.to_dict(SeqIO.parse(handle, "gb"))


async def get_qualifiers(seq) -> dict:
    """
    Get relevant qualifiers in a Genbank record

    :param seq: SeqIO features object for a particular accession
    :return: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    """
    f = [feature for feature in seq if feature.type == "source"]
    isolate_data = {}

    for feature in f:
        for qualifier in feature.qualifiers:
            isolate_data[qualifier] = feature.qualifiers.get(qualifier)
    return isolate_data


async def find_isolate(isolate_data) -> str:
    """
    Determine the source type in a Genbank record

    :param isolate_data: Dictionary containing qualifiers in a features section of a Genbank record
    :return:
    """
    isolate_type = None
    for qualifier in ["isolate", "strain"]:
        if qualifier in isolate_data:
            isolate_type = qualifier
            break

    return isolate_type


async def store_sequence(path, accession, data, unique_ids) -> Union[str, None]:
    """
    Creates a new sequence file for a given isolate

    :param path: Path to a isolate folder
    :param accession: Genbank record object for a given accession
    :param data: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    :param unique_ids: Set containing unique Virtool ids for all sequences in a reference
    :return: The newly generated unique id
    """
    sequences = await get_sequences(path)

    # check if accession doesn't already exist
    if accession.id in sequences:
        return

    # generate new sequence id
    new_id = random_alphanumeric(8, False, unique_ids)

    seq_file = {"_id": new_id,
                "accession": accession.id,
                "definition": accession.description,
                "host": data.get("host")[0] if data.get("host") is not None else None,
                "sequence": str(accession.seq)
                }

    async with aiofiles.open(os.path.join(path, new_id + ".json"), "w") as f:
        await f.write(json.dumps(seq_file, indent=4))
        await f.seek(0)

    return new_id


async def store_isolate(path, source_name, source_type, unique_ids) -> str:
    """
    Creates a new isolate folder for an OTU

    :param path: Path to an OTU
    :param source_name: Assigned source name for an accession
    :param source_type: Assigned source type for an accession
    :param unique_ids: Set containing unique Virtool ids for all isolates in a reference
    :return: A newly generated unique id
    """

    new_id = random_alphanumeric(8, False, unique_ids)
    os.mkdir(os.path.join(path, new_id))

    new_isolate = {
        "id": new_id,
        "source_type": source_type,
        "source_name": source_name,
        "default": False
    }

    async with aiofiles.open(os.path.join(path, new_id, "isolate.json"), "w") as f:
        await f.write(json.dumps(new_isolate, indent=4))

    return new_id


def get_cache() -> dict:
    """
    Fetches the local cache containing taxon ids mapped to their accessions

    :return: A dictionary that maps taxon ids to found accessions
    """
    if not os.path.isdir(".cli"):
        return {}

    with open(".cli/accessions_cache.json", "r") as f:
        cache = json.load(f)

        return cache


def write_cache(accessions):
    """
    Cache a mapping of taxon ids to all accessions found.

    :param accessions: Dictionary containing an updated mapping of taxids to their accessions
    """
    if not os.path.isdir(".cli"):
        os.mkdir(".cli")

    with open(".cli/accessions_cache.json", "w") as f:
        f.seek(0)
        json.dump(accessions, f, indent=4)


def random_alphanumeric(length: int = 6, mixed_case: bool = False, excluded: Union[None, Iterable[str]] = None) -> str:
    """
    Generates a random string composed of letters and numbers.

    :param length: the length of the string.
    :param mixed_case: included alpha characters will be mixed case instead of lowercase
    :param excluded: strings that may not be returned.
    :return: a random alphanumeric string.
    """
    excluded = set(excluded or list())

    characters = digits + (ascii_letters if mixed_case else ascii_lowercase)

    candidate = "".join([choice(characters) for _ in range(length)])

    if candidate not in excluded:
        return candidate

    return random_alphanumeric(length=length, excluded=excluded)


def run(src):
    """
    Runs the asynchronous routines to find new isolates for all OTU in a reference

    :param src: Path to a reference directory
    """
    asyncio.run(isolate(src))
