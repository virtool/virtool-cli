import asyncio
import json
import pathlib
from concurrent.futures.thread import ThreadPoolExecutor
from typing import Optional

import aiofiles
import aiojobs
from Bio import Entrez
from rich.console import Console

from virtool_cli.utils import get_otu_paths, NCBI_REQUEST_INTERVAL


async def taxid(src_path: pathlib.Path, force_update: bool):
    """
    Asynchronously finds taxon ids for all OTU in a given src directory and
    writes them to each respective otu.json file.

    :param src_path: Path to a given reference directory
    :param force_update: Flag to force update all of the OTU's taxids in a reference
    """
    scheduler = await aiojobs.create_scheduler()
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor())

    paths = get_otu_paths(src_path)
    coros = []
    otu_paths = {}

    for path in paths:
        name = await get_name_from_path(path, force_update)
        if name is not None:
            otu_paths[name] = path
            coros.append(name)

    console = Console()
    results = list()

    q = asyncio.Queue()

    # Put coroutines into the scheduler as long as its active number of jobs doesn't exceed the concurrency limit
    while len(coros) != 0:
        if scheduler.active_count < scheduler.limit:
            name = coros.pop()

            await scheduler.spawn(fetch_taxid_call(name, q, console))

            await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    # Pulling results from queue. Don't stop checking until the queue is empty and the scheduler has no active jobs.
    while True:
        try:
            results.append(q.get_nowait())
        except asyncio.QueueEmpty:
            if not scheduler.active_count:
                break

            await asyncio.sleep(0.01)

        await scheduler.close()

    names = []
    taxids = []

    for name, taxid in results:
        names.append(name)
        if taxid:
            taxids.append(taxid)
        update_otu(taxid, otu_paths[name])

    console.print(
        f"\nRetrieved {len(taxids)} taxids for {len(names)} OTUs", style="green"
    )


def fetch_taxid(name: str) -> int:
    """
    Searches the NCBI taxonomy database for a given OTU name and fetches and returns
    its taxon id. If a taxon id was not found the function should return None.

    :param name: Name of a given OTU
    :return: Taxonomy id for the given OTU
    """
    handle = Entrez.esearch(db="taxonomy", term=name)
    record = Entrez.read(handle)

    try:
        taxid = int(record["IdList"][0])
    except IndexError:
        taxid = None

    return taxid


async def fetch_taxid_call(name: str, q: asyncio.queues.Queue, console: Console):
    """
    Handles calling asynchronous taxon id retrievals and updating task progress.
    Puts results in a asyncio Queue.

    :param name: Name of a given OTU
    :param q: asyncio Queue used to put results into
    :param console: Rich console object used for logging
    """
    taxid = await asyncio.get_event_loop().run_in_executor(None, fetch_taxid, name)

    await log_results(name, taxid, console)

    await q.put((name, taxid))


async def log_results(name: str, taxid: int, console: Console):
    """
    Logs results of the taxid fetch from the NCBI taxonomy database

    :param name: Name of a given OTU
    :param taxid: Taxid for a given OTU if found, else None
    :param console: Rich console object used for logging
    """
    if taxid:
        console.print(f"[green]✔ {name}")
        console.print(f"[green]  Found taxid {taxid}")
    else:
        console.print(f"[red]✘ {name}")
        console.print("[red]  Could not find taxid")

    console.print()


def update_otu(taxid: int, path: pathlib.Path):
    """
    Updates a otu.json's taxid key with either a taxon id or None
    depending on whether an id was able to be retrieved.

    :param taxid: Taxid for a given OTU if found, else None
    :param path: Path to a given OTU in the reference directory
    """
    with open(path / "otu.json", "r+") as f:
        otu = json.load(f)
        otu["taxid"] = taxid if taxid else taxid
        f.seek(0)
        json.dump(otu, f, indent=4)


async def get_name_from_path(path: pathlib.Path, force_update: bool) -> Optional[str]:
    """
    Given a path to an OTU, returns the name of the OTU. If a taxon id already exists
    for the OTU then it returns None.

    :param path: Path to a given OTU
    :param force_update: Boolean flag used to determine if the function should return the name of the OTU even if it
    already has a taxid
    :return: The OTU's name if it doesn't have a taxid assigned or force_update is True, else None
    """
    async with aiofiles.open(path / "otu.json", "r") as f:
        otu = json.loads(await f.read())

        try:
            if otu["taxid"] is not None and not force_update:
                return None
            else:
                return otu["name"]
        except KeyError:
            return otu["name"]


def run(src_path: pathlib.Path, force_update: bool):
    """
    Creates and runs the event loop to asynchronously find taxon ids for all OTU in a src directory

    :param src_path: Path to a given reference directory
    :param force_update: Flag to force update all of the OTU's taxids in a reference
    :return:
    """
    asyncio.run(taxid(src_path, force_update))
