import asyncio
import json
import os
from concurrent.futures.thread import ThreadPoolExecutor
from random import choice
from string import ascii_letters, ascii_lowercase, digits
from typing import Union, Iterable
from urllib.error import HTTPError

import aiofiles
import aiojobs
from Bio import Entrez, SeqIO
from rich.console import Console

from virtool_cli.utils import get_otu_paths, get_taxid_map, get_isolates, get_sequences, get_unique_ids

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

REQUEST_INTERVAL = 0.3 if Entrez.email and Entrez.api_key else 0.8


async def isolate(src):
    scheduler = await aiojobs.create_scheduler(limit=10)
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor())

    q = asyncio.Queue()

    existing_accessions = check_accessions_cache()

    # get paths for all OTU in the directory
    paths = get_otu_paths(src)

    # get mapping of all OTU paths to their taxid
    taxid_otu_path_map = await get_taxid_map(paths)

    unique_ids = await get_unique_ids(paths)

    for path in paths:
        # only fetch OTU that have a taxid
        taxid = str(taxid_otu_path_map[path])

        if taxid is None:
            continue

        accessions = existing_accessions.get(taxid)

        await scheduler.spawn(fetch_otu_isolates(taxid, path, accessions, unique_ids, q))

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

    new_cache = dict()
    for taxid, accessions in results_to_cache:
        new_cache[taxid] = accessions

    existing_accessions.update(new_cache)

    cache_accessions(existing_accessions)


async def fetch_otu_isolates(taxid, path, accessions, unique_ids, queue):
    # get existing isolates
    isolates = await get_isolates(path)

    records, new_accessions = await asyncio.get_event_loop().run_in_executor(None, get_records, accessions, taxid)

    if records is not None:
        for accession in [accession for accession in records.values() if accession.seq]:
            isolate_data = await get_qualifiers(accession.features)

            # try to find isolate type and name automatically
            isolate_type = None
            for qualifier in ["isolate", "strain"]:
                if qualifier in isolate_data:
                    isolate_type = qualifier
                    break

            if not isolate_type:
                continue

            # see if isolate directory already exists and create it if needed
            isolate_name = isolate_data.get(isolate_type)[0]

            if isolate_name not in isolates:
                new_id = await store_isolate(path, isolate_name, isolate_type, unique_ids)
                unique_ids.add(new_id)
                isolates[isolate_name] = new_id

            isolate_path = os.path.join(path, isolates.get(isolate_name))
            new_id = await store_sequence(isolate_path, accession, isolate_data, unique_ids)
            unique_ids.add(new_id)

            await queue.put((taxid, new_accessions))


def get_records(accessions, taxid):
    console = Console()

    if accessions is None:
        accessions = []
        record = Entrez.read(Entrez.elink(
            dbfrom="taxonomy", db="nucleotide", id=taxid, idtype="acc"))

        for linksetdb in record[0]["LinkSetDb"][0]["Link"]:
            accessions.append(linksetdb["Id"])

    try:
        records = fetch_records(accessions)
        console.print(f"Found isolate data for {taxid}", style="green")
    except HTTPError:
        console.print(f"Could not find isolate data for {taxid}", style="red")
        return None, None

    return records, accessions


def fetch_records(accessions):
    handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb",
                           retmode="text")
    return SeqIO.to_dict(SeqIO.parse(handle, "gb"))


async def get_qualifiers(seq):
    f = [feature for feature in seq if feature.type == "source"]
    isolate_data = {}

    for feature in f:
        for qualifier in feature.qualifiers:
            isolate_data[qualifier] = feature.qualifiers.get(qualifier)
    return isolate_data


def cache_accessions(accessions):
    """Cache a mapping of taxon ids to all accessions found."""
    if not os.path.isdir(".cli"):
        os.mkdir(".cli")

    with open(".cli/accessions_cache.json", "w") as f:
        f.seek(0)
        json.dump(accessions, f, indent=4)


async def store_sequence(path, accession, data, unique_ids):
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


async def store_isolate(path, source_name, source_type, unique_ids):
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


def check_accessions_cache():
    if not os.path.isdir(".cli"):
        return {}

    with open(".cli/accessions_cache.json", "r") as f:
        cache = json.load(f)

        return cache


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
    asyncio.run(isolate(src))
