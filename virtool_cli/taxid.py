import os
import time
import json
import asyncio
from rich.console import Console
from rich.progress import track
from rich.progress import Progress
from Bio import Entrez
from concurrent.futures.thread import ThreadPoolExecutor

Entrez.email = os.environ["NCBI_EMAIL"]
API_KEY = os.environ["NCBI_API_KEY"]
progress = Progress()


async def run(src_path: str):
    executor = ThreadPoolExecutor(max_workers=10)
    loop = asyncio.get_event_loop()

    paths = get_paths(src_path)

    coros = list()
    console = Console()
    otus = {}

    for path in track(paths[:100], description="[green]Retrieving..."):
        name = get_name_from_path(path)

        if name:
            coro = await fetch_taxid_call(executor, name, path, console, loop)
            coros.append(coro)
            otus[name] = path
            await asyncio.sleep(0.2)

    results = await asyncio.gather(*coros)

    names = []
    taxids = []

    for name, taxid in results:
        names.append(name)
        if taxid:
            taxids.append(taxid)

        update_otu(name, taxid, otus[name])
        await asyncio.sleep(0.05)

    console.print(
        f"Retrieved {len(taxids)} taxids for {len(names)} OTUs", style="green")


def fetch_taxid(name, path):
    handle = Entrez.esearch(db="taxonomy", term=name, api_key=API_KEY)
    record = Entrez.read(handle)

    try:
        taxid = record["IdList"][0]
    except IndexError:
        taxid = None

    return name, taxid


async def fetch_taxid_call(executor, name, path, console, loop):
    console.print(f"    {name}")
    coro = loop.run_in_executor(executor, fetch_taxid, name, path)

    return coro


def update_otu(name, taxid, path):
    console = Console()
    with open(os.path.join(path, "otu.json"), 'r+') as f:
        otu = json.load(f)

        if taxid:
            console.print(f"Retrieved {name} {taxid} :heavy_check_mark:", style="green")
        else:
            console.print(f"Could not receive {name} :x:", style="red")
        otu["taxid"] = taxid
        f.seek(0)
        json.dump(otu, f, indent=4)


def get_paths(src_path: str):
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


def get_name_from_path(path):
    with open(os.path.join(path, "otu.json"), 'r') as f:
        otu = json.load(f)

        try:
            if otu["taxid"]:
                return None
        except KeyError:
            return otu["name"]


def taxid(src_path):
    loop = asyncio.get_event_loop()
    loop.run_until_complete(run(src_path))
