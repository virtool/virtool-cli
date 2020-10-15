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

from virtool_cli.utils import get_paths, get_taxids, get_isolates

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

qualifiers = ["isolate", "strain"]

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")
bad = {"bad_requests": []}


async def isolate(src):
    scheduler = await aiojobs.create_scheduler()
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor())

    global accessions
    accessions = check_accessions_cache() if check_accessions_cache() is not None else {}

    # get paths for all OTU in the directory
    paths = get_paths(src)

    # get mapping of all OTU paths to their taxid
    taxids = await get_taxids(paths)

    for path in paths[:10]:
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
    
    with open("bad_requests.json", "w") as f:
        f.seek(0)
        json.dump(bad, f, indent=4)


async def run(taxid, path):
    global accessions

    # get mapping of all isolates to their folder name
    isolates = await get_isolates(path)

    records = await asyncio.get_event_loop().run_in_executor(None, fetch, accessions, taxid)

    if records is not None:
        no_data = 0
        for seq in [seq for seq in records.values() if seq.seq]:
            isolate_data = await get_qualifiers(seq.features)
            if isolate_data["isolate"] is None and isolate_data["strain"] is None:
                no_data += 1
        print(
            f"Isolate source data not found in {no_data} out of {len(records)} records")
    else:
        bad["bad_requests"].append(taxid)


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
        for qualifier in qualifiers:
            isolate_data[qualifier] = feature.qualifiers.get(qualifier)
    return isolate_data


def cache_accessions(accessions):
    if not os.path.isdir(".cli"):
        os.mkdir(".cli")

    with open(".cli/accessions_cache.json", "w") as f:
        f.seek(0)
        json.dump(accessions, f, indent=4)


def cache_records(taxid, records):
    with open(".cli/accessions_cache.json", "r+") as f:
        cache = json.load(f)

        # can't cache features since they are SeqRecord objects
        for record in records.values():
            new_record = {"definition": record.description,
                          "url": f"https://www.ncbi.nlm.nih.gov/nuccore/{record.id}", "sequence": str(record.seq),
                          "features": [str(feature) for feature in record.features if feature.type == "source"]}
            cache[taxid][record.id] = new_record
        f.seek(0)
        json.dump(cache, f, indent=4)


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
    asyncio.run(isolate("tests/files/src"))
