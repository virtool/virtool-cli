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

import virtool_cli.taxid

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
    paths = virtool_cli.taxid.get_paths(src)

    for path in paths:
        # only fetch OTU that have a taxid
        taxid = await get_taxid(path)
        if not taxid:
            continue

        # determine retmax
        handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[orgn]", retmax=0)
        record = Entrez.read(handle)
        search_max = record["Count"]

        handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[orgn]", retmax=search_max)
        record = Entrez.read(handle)
        alist = record["IdList"]

        # approach 1
        handle = Entrez.efetch(db="nuccore", id=alist, rettype="gb", retmode="text")
        record = SeqIO.parse(handle, "gb")

        for seq in record:
            print()
            for feature in seq.features:
                if feature.type == "source":
                    print(feature.qualifiers['isolate'], seq.seq)


async def run():
    pass


async def get_taxid(path):
    async with aiofiles.open(os.path.join(path, "otu.json"), "r") as f:
        otu = json.loads(await f.read())
        return otu["taxid"] if "taxid" in otu else None


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
