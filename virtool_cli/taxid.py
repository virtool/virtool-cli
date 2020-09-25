"""
TODO:
- add docstrings
- look at aiojobs docs
- look at docs about asyncio.Queue
bonus: figure out why it hangs at the very end
"""

import asyncio
import json
import os
from concurrent.futures.thread import ThreadPoolExecutor

import aiofiles
import aiojobs
from Bio import Entrez
from rich.console import Console
from rich.progress import Progress

Entrez.email = os.environ["NCBI_EMAIL"]
API_KEY = os.environ["NCBI_API_KEY"]


async def run(src_path: str, force_update: bool):
    scheduler = await aiojobs.create_scheduler(limit=5)
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor(max_workers=5))

    paths = get_paths(src_path)

    console = Console()
    results = list()

    q = asyncio.Queue()
    otu_paths = {}

    with Progress("[progress.description]{task.description}", "{task.fields[result]}") as progress:
        # Submit jobs to scheduler and create progress tasks.
        for path in paths:
            name = await get_name_from_path(path, force_update)

            if name:
                otu_paths[name] = path
                task = progress.add_task(description=f"Retrieving taxon ID for {name}", result="")

                await scheduler.spawn(fetch_taxid_call(name, progress, q, task))
                await asyncio.sleep(0.2)

        # Pulling results from queue. Don't stop checking until the queue is empty and the scheduler has no active jobs.
        while True:
            try:
                results.append(q.get_nowait())
            except asyncio.QueueEmpty:
                if not scheduler.active_count:
                    break

            await asyncio.sleep(0.3)

        await scheduler.close()

    names = []
    taxids = []

    for name, taxid in results:
        names.append(name)
        if taxid:
            taxids.append(taxid)
        update_otu(taxid, otu_paths[name])

    console.print(f"\nRetrieved {len(taxids)} taxids for {len(names)} OTUs", style="green")


def fetch_taxid(name):
    handle = Entrez.esearch(db="taxonomy", term=name, api_key=API_KEY)
    record = Entrez.read(handle)

    try:
        taxid = record["IdList"][0]
    except IndexError:
        taxid = None

    return name, taxid


async def fetch_taxid_call(name, progress, q, task):
    name, taxid = await asyncio.get_event_loop().run_in_executor(None, fetch_taxid, name)

    if taxid is None:
        description = f"[red]Retrieving taxon ID for {name}"
        result = f"[red]:x: Not found"
    else:
        description = f"[green]Retrieving taxon ID for {name}"
        result = f"[green]:heavy_check_mark: {taxid}"

    progress.update(task, description=description, result=result)

    await q.put((name, taxid))


def update_otu(taxid, path):
    with open(os.path.join(path, "otu.json"), 'r+') as f:
        otu = json.load(f)
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


# Changed this to use aiofiles
async def get_name_from_path(path, force_update):
    async with aiofiles.open(os.path.join(path, "otu.json"), "r") as f:
        otu = json.loads(await f.read())

        try:
            if otu["taxid"] is not None and not force_update:
                return None
            else:
                return otu["name"]
        except KeyError:
            return otu["name"]


def taxid(src_path, force_update):
    asyncio.run(run(src_path, force_update))
