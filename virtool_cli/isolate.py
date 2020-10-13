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

qualifiers = ["isolate", "strain"]

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")


async def isolate(src):
    scheduler = await aiojobs.create_scheduler(limit=5)
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor(max_workers=5))

    # get paths for all OTU in the directory
    paths = get_paths(src)

    # get mapping of all OTU paths to their taxid
    taxids = await get_taxids(paths)
    unknown = 0
    for path in paths[:200]:

        # only fetch OTU that have a taxid
        if taxids[path] is None:
            continue

        # get mapping of all isolates to their folder name
        isolates = await get_isolates(path)

        accessions = await check_accessions_cache(taxids[path])

        if not accessions:
            accessions = {}
            record = Entrez.read(Entrez.elink(
                dbfrom="taxonomy", db="nucleotide", id=taxids[path], idtype="acc"))

            for linksetdb in record[0]["LinkSetDb"][0]["Link"]:
                accessions[linksetdb["Id"]] = None

            await cache_accessions(str(taxids[path]), accessions)
        try:
            handle = Entrez.efetch(db="nucleotide", id=[accession for accession in accessions.keys()], rettype="gb",
                                   retmode="text")
        except HTTPError:
            unknown += 1

        # this is somehow exponentially slower than the fetch itself
        records = SeqIO.to_dict(SeqIO.parse(handle, "gb"))

        no_data = 0
        for seq in [seq.features for seq in records.values()]:
            isolate_data = await get_qualifiers(seq)
            if isolate_data["isolate"] is None and isolate_data["strain"] is None:
                no_data += 1
        print(
            f"Isolate source data not found in {no_data} out of {len(records)} records")

    print(f"{unknown} taxids were not found in the nucleotide database")


async def run():
    pass


async def get_qualifiers(seq):
    f = [feature for feature in seq if feature.type == "source"]
    isolate_data = {}

    for feature in f:
        for qualifier in qualifiers:
            isolate_data[qualifier] = feature.qualifiers.get(qualifier)
    return isolate_data


async def cache_accessions(taxid, accessions):
    if not os.path.isdir(".cli"):
        os.mkdir(".cli")
        cache = {taxid: accessions}
    else:
        async with aiofiles.open(".cli/accessions_cache.json", "r") as f:
            cache = json.loads(await f.read())
            cache[taxid] = accessions

    async with aiofiles.open(".cli/accessions_cache.json", "w") as f:
        f.seek(0)
        await f.write(json.dumps(cache, indent=4))


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


async def check_accessions_cache(taxid):
    if not os.path.isdir(".cli"):
        return None

    async with aiofiles.open(".cli/accessions_cache.json", "r") as f:
        cache = json.loads(await f.read())

        return cache.get(str(taxid))


async def check_records_cache():
    pass


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
    if os.path.exists(".cli"):
        shutil.rmtree(".cli")
    asyncio.run(isolate("tests/files/src"))
