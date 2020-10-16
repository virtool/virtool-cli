import asyncio
import json
import os
from concurrent.futures.thread import ThreadPoolExecutor
from random import choice

import shutil
from string import ascii_letters, ascii_lowercase, digits
from typing import Union, Iterable
from urllib.error import HTTPError
import aiofiles
import aiojobs
from Bio import Entrez, SeqIO

from virtool_cli.utils import *

ISOLATE_KEYS = [
    "id",
    "source_type",
    "source_name",
    "default"
]

SEQUENCE_KEYS = [
    "_id",
    "accession",
    "definition",
    "host",
    "sequence"
]

accessions = {}
unique_ids = {}

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")


async def isolate(src):
    scheduler = await aiojobs.create_scheduler()
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor())

    global accessions
    accessions = check_accessions_cache() if check_accessions_cache() is not None else {}

    # get paths for all OTU in the directory
    paths = get_paths(src)

    # get mapping of all OTU paths to their taxid
    taxids = await get_taxids(paths)

    global unique_ids
    unique_ids = await get_unique_ids(paths)

    for path in paths:
        # only fetch OTU that have a taxid
        if taxids[path] is None:
            continue

        await scheduler.spawn(run(str(taxids[path]), path))
        if Entrez.email is None or Entrez.api_key is None:
            await asyncio.sleep(0.8)
        else:
            await asyncio.sleep(0.2)

    while True:
        if not scheduler.active_count:
            await scheduler.close()
            break
        await asyncio.sleep(0.1)

    cache_accessions(accessions)


async def run(taxid, path):

    isolates = await get_isolates(path)

    records = await asyncio.get_event_loop().run_in_executor(None, fetch, accessions, taxid)
    
    if records is not None:
        for accession in [accession for accession in records.values() if accession.seq]:
            isolate_data = await get_qualifiers(accession.features)

            # try to find isolate type and name automatically
            found = None
            for qualifier in ["isolate", "strain"]:
                if qualifier in isolate_data:
                    found = qualifier
                    break
            if found:
                # see if isolate directory already exists
                if isolate_data.get(found)[0] in isolates:
                    isolate_name = isolate_data.get(found)[0]
                    isolate_path = os.path.join(path, isolates.get(isolate_name))
                    await store_sequence(isolate_path, accession, isolate_data)
                else:
                    pass
            # else pass manual input to user
            else:
                pass


def fetch(accessions, taxid):
    if not accessions.get(taxid):
        record = Entrez.read(Entrez.elink(
            dbfrom="taxonomy", db="nucleotide", id=taxid, idtype="acc"))
        acc_dict = {}

        for linksetdb in record[0]["LinkSetDb"][0]["Link"]:
            acc_dict[linksetdb["Id"]] = None
        accessions[taxid] = acc_dict
    try:
        handle = Entrez.efetch(db="nucleotide", id=[accession for accession in accessions[taxid]], rettype="gb",
                               retmode="text")
    except HTTPError:
        return None

    return SeqIO.to_dict(SeqIO.parse(handle, "gb"))


async def get_qualifiers(seq):
    f = [feature for feature in seq if feature.type == "source"]
    isolate_data = {}

    for feature in f:
        for qualifier in feature.qualifiers:
            isolate_data[qualifier] = feature.qualifiers.get(qualifier)
    return isolate_data


def cache_accessions(accessions):
    """Cache a mapping of taxon ids to all accessions found"""
    if not os.path.isdir(".cli"):
        os.mkdir(".cli")

    with open(".cli/accessions_cache.json", "w") as f:
        f.seek(0)
        json.dump(accessions, f, indent=4)


async def store_sequence(path, new_seq, data):

    sequences = await get_sequences(path)

    # check if accession doesn't already exist
    if new_seq.id in sequences:
        return

    # generate new sequence id
    new_id = random_alphanumeric(8, False, unique_ids)

    await update_ids(new_id)

    seq_file = {"_id": new_id,
                "accession": new_seq.id,
                "definition": new_seq.description,
                "host": data.get("host")[0] if data.get("host") is not None else None,
                "sequence": str(new_seq.seq)
                }

    async with aiofiles.open(os.path.join(path, new_id + ".json"), "w") as f:
        await f.write(json.dumps(seq_file, indent=4))


async def update_ids(id):
    global unique_ids

    unique_ids.add(id)


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


if __name__ == '__main__':
    asyncio.run(isolate("tests/files/src_d"))
