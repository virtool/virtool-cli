import asyncio
import json
import os
from concurrent.futures.thread import ThreadPoolExecutor
from random import choice
from string import ascii_letters, ascii_lowercase, digits
from typing import Union, Iterable

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

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")


async def isolate(src):
    scheduler = await aiojobs.create_scheduler(limit=5)
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor(max_workers=5))

    # get paths for all OTU in the directory
    paths = get_paths(src)
    taxids = await get_taxids(paths)
    cache = {}
    for path in taxids.keys():
        # only fetch OTU that have a taxid
        if taxids[path] is None:
            continue

        # get mapping of all isolates to their folder name
        isolates = await get_isolates(path)

        # approach 1
        accessions = await check_cache(taxids[path])

        if not accessions:
            accessions = []
            record = Entrez.read(Entrez.elink(dbfrom="taxonomy", db="nucleotide", id=taxids[path]))
            for linksetdb in record[0]["LinkSetDb"][0]["Link"]:
                accessions.append(linksetdb["Id"])

            cache_taxid(taxids[path], accessions)

        handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="text")
        record = SeqIO.parse(handle, "gb")

        count = []
        for seq in record:
            for feature in seq.features:
                if feature.type == "source":
                    if "isolate" not in feature.qualifiers:
                        count.append(seq)

        actual = []
        print("Accessions found:", len(accessions))
        print("Accessions with no isolate in feature qualifier:", len(count))
        for seq in count:
            if "isolate" not in seq.description:
                actual.append(seq)
        print("Accessions with no isolate in feature qualifier AND no isolate in definition:", len(actual))


async def run():
    pass


def cache_taxid(taxid, accessions):
    if not os.path.isdir(".cli"):
        os.mkdir(".cli")
        with open(".cli/accessions_cache.json", "w") as f:
            json.dump({}, f, indent=4)

    with open(".cli/accessions_cache.json", "r+") as f:
        cache = json.load(f)
        cache[taxid] = accessions
        f.seek(0)
        json.dump(cache, f, indent=4)


async def check_cache(taxid):
    if not os.path.isdir(".cli"):
        return None

    async with aiofiles.open(".cli/accessions_cache.json", "r") as f:
        cache = json.loads(await f.read())
        try:
            return cache[str(taxid)]
        except KeyError:
            return None


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
    asyncio.run(isolate("tests/files/src_b"))
