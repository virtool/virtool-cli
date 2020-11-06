import asyncio
import json
import os
from concurrent.futures.thread import ThreadPoolExecutor
from random import choice
from string import ascii_letters, ascii_lowercase, digits
from typing import Union, Iterable, Tuple, Optional
from urllib.error import HTTPError

import aiofiles
import aiojobs
from Bio import Entrez, SeqIO
from rich.console import Console

from virtool_cli.utils import get_otu_paths, get_otus, get_isolates, get_sequences, \
    get_unique_ids, NCBI_REQUEST_INTERVAL


async def isolate(src):
    """
    Runs routines to find new isolates for OTU in a reference directory and writes newfound accessions to local cache

    :param src: Path to a given reference directory
    """
    scheduler = await aiojobs.create_scheduler(limit=10)
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor())

    console = Console()

    q = asyncio.Queue()

    existing_accessions = get_cache()

    # get paths for all OTU in the directory
    paths = get_otu_paths(src)

    # get mapping of all OTU paths to its OTU JSON file
    otu_path_map = await get_otus(paths)

    for path in paths:
        taxid = otu_path_map[path].get("taxid")
        name = otu_path_map[path].get("name")

        if taxid is None:
            console.print(f"✘ [red]{name}\n"
                          f"  [red]No taxid assigned to OTU\n")
            continue

        accessions = existing_accessions.get(str(taxid))
        try:
            await scheduler.spawn(fetch_otu_isolates(str(taxid), name, path, accessions, paths, q))
        except KeyboardInterrupt:
            await scheduler.close()
            break

        await asyncio.sleep(NCBI_REQUEST_INTERVAL)

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


async def fetch_otu_isolates(taxid, name, path, accessions, paths, queue):
    """
    Fetch new isolates for a given OTU, creating a directory for each

    :param taxid: The taxon id for a given OTU
    :param name: Name of a given OTU
    :param path: Path to a given OTU in a reference directory
    :param accessions: List of accessions for the given taxon id from local cache
    :param paths: List of all paths to every OTU in a reference
    :param queue: Asyncio queue to put newfound accessions onto for caching purposes
    """

    isolates = await get_isolates(path)

    records, new_accessions = await asyncio.get_event_loop().run_in_executor(None, get_records, accessions, taxid)

    isolate_ids, sequence_ids = await get_unique_ids(paths)

    console = Console()

    if records is None:
        console.print(f"[red]✘ {name} ({taxid})")
        console.print("  [red]Found 0 isolates, could not link taxid to nucleotide records\n")
        return

    new_isolates = dict()

    for accession in [accession for accession in records.values() if accession.seq]:
        accession_data = await get_qualifiers(accession.features)

        # try to find isolate type and name automatically
        isolate_type = await find_isolate(accession_data)

        if isolate_type is None:
            continue

        isolate_name = accession_data.get(isolate_type)[0]

        # see if isolate directory already exists and create it if needed
        if isolate_name not in isolates:
            new_isolates[isolate_name] = isolate_type

            new_id = await store_isolate(path, isolate_name, isolate_type, isolate_ids)
            isolate_ids.add(new_id)
            isolates[isolate_name] = new_id

        # create new sequence file
        isolate_path = os.path.join(path, isolates.get(isolate_name))
        new_id = await store_sequence(isolate_path, accession, accession_data, sequence_ids)
        sequence_ids.add(new_id)

        await queue.put((taxid, new_accessions))

    await log_results(name, taxid, new_isolates, console)


def get_records(accessions, taxid) -> Tuple[Optional[dict], list]:
    """
    Uses an OTU's taxid to find accessions numbers if needed, then fetches Genbank records for each accession number

    :param accessions: A list of accession ids to retrieve or None to retrieve accession numbers for the given taxid
    :param taxid: The taxid to find accessions for
    :return: A dictionary with all the accession records and a list containing all accession ids for caching purposes
    """

    if accessions is None:
        accessions = get_accession_numbers(taxid)

    records = fetch_records(accessions)

    return records, accessions


def get_accession_numbers(taxid) -> list:
    """
    Links an OTU's taxid to all associated accession numbers in the Nucleotide database

    :param taxid: The taxid to find accessions for
    :return: List containing accession numbers to fetch records for
    """
    accessions = []

    record = Entrez.read(Entrez.elink(dbfrom="taxonomy", db="nucleotide",
                                      id=taxid, idtype="acc"))

    for linksetdb in record[0]["LinkSetDb"][0]["Link"]:
        accessions.append(linksetdb["Id"])

    return accessions


def fetch_records(accessions) -> Optional[dict]:
    """
    Fetches a Genbank record for each accession number in a given list, then parses all of the records into a
    SeqIO dictionary that can be used to access specific fields in each record

    :param accessions: List of accession numbers to be fetched
    :return: A SeqIO dictionary object
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=accessions,
                               rettype="gb", retmode="text")
    except HTTPError:
        return None

    return SeqIO.to_dict(SeqIO.parse(handle, "gb"))


async def log_results(name, taxid, new_isolates, console):
    """
    Log isolate discovery results to console

    :param name: Name of the OTU to log
    :param taxid: Taxid of the OTU to log
    :param new_isolates: Dictionary containing newly found isolates mapped to their type
    :param console: Rich console object used to log results
    """
    new_isolate_count = len(new_isolates)

    if new_isolate_count == 0:
        console.print(f"[red]✘ {name} ({taxid})\n  Found {new_isolate_count} new isolates\n")
    else:
        console.print(f"[green]✔ {name} ({taxid})\n  Found {new_isolate_count} new isolates:")
        for isolate_name, isolate_type in new_isolates.items():
            console.print(f"[green]    - {isolate_type} {isolate_name}")
        console.print()


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

    # strip accession version
    accession_id = accession.id.split(".")[0]

    # check if accession doesn't already exist
    if accession_id in sequences:
        return

    # generate new sequence id
    new_id = random_alphanumeric(8, False, unique_ids)

    seq_file = {"_id": new_id,
                "accession": accession_id,
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


def update_cache(existing_cache, results) -> dict:
    """
    Updates the existing accessions cache with fetch results

    :param existing_cache: The existing cache from the .cli folder
    :param results: List of tuples pulled out of the asyncio Queue containing each taxid and their accessions
    """
    cache_update = dict(results)

    return {
        **existing_cache,
        **cache_update
    }


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
