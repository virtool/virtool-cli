import asyncio
import json
import os
from typing import Union
from concurrent.futures.thread import ThreadPoolExecutor

import aiofiles
import aiojobs
from Bio import Entrez
from virtool_cli.utils import get_otu_paths, NCBI_REQUEST_INTERVAL
from rich.console import Console
from rich.progress import BarColumn, Progress, TimeRemainingColumn

missed_otus = {"otus": []}


async def taxid(src_path: str, force_update: bool):
    """
    Asynchronously finds taxon ids for all OTU in a given src directory and
    writes them to each respective otu.json file.

    Parameters:
        src_path (str): Path to database src directory
        force_update (bool): Flag to indicate whether all OTUs should be searched for regardless
                             if they already have a taxon id


    """
    scheduler = await aiojobs.create_scheduler(limit=5)
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor(max_workers=5))

    paths = get_otu_paths(src_path)
    checked_names = []
    otu_paths = {}

    for path in paths:
        name = await get_name_from_path(path, force_update)
        if name is not None:
            otu_paths[name] = path
            checked_names.append(name)

    console = Console()
    results = list()

    q = asyncio.Queue()

    progress = Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        "[magenta]{task.completed} of {task.total} ids searched",
        TimeRemainingColumn()
    )

    with progress:
        task = progress.add_task("Retrieving...", total=len(checked_names))
        for name in checked_names:
            progress.update(task, description=f"Retrieving {name}")

            await scheduler.spawn(fetch_taxid_call(name, progress, q, task))

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

    console.print(f"\nRetrieved {len(taxids)} taxids for {len(names)} OTUs", style="green")
    if len(missed_otus["otus"]):
        console.print(f"OTUs that could not be retrieved have been logged to missed_otus.json", style="green")
        with open("missed_otus.json", 'w') as f:
            json.dump(missed_otus, f, indent=4)


def fetch_taxid(name: str) -> (str, str):
    """
    Searches the NCBI taxonomy database for a given OTU name and fetches and returns
    its taxon id. If a taxon id was not found the function should return None.

    Parameters:
        name (str): Name of a OTU

    Returns:
        Taxon id for a given OTU
    """
    handle = Entrez.esearch(db="taxonomy", term=name)
    record = Entrez.read(handle)

    try:
        taxid = record["IdList"][0]
    except IndexError:
        taxid = None
        missed_otus["otus"].append(name)

    return taxid


async def fetch_taxid_call(name: str, progress: Progress, q: asyncio.queues.Queue, task: int):
    """
    Handles calling asynchronous taxon id retrievals and updating task progress.
    Puts results in a asyncio Queue.

    Parameters:
        name (str): name of a OTU
        progress (Progress): Progress object for current tasks that can be updated
        q (asyncio.queues.Queue): Queue containing current results of tasks
        task (int): ID for a given progress task

    """
    taxid = await asyncio.get_event_loop().run_in_executor(None, fetch_taxid, name)

    if taxid is None:
        description = f"[red]Could not retrieve taxon ID for {name}"
        result = f"[red]:x: Not found"
    else:
        description = f"[green]Retrieved taxon ID for {name}"
        result = f"[green]:heavy_check_mark: {taxid}"

    progress.print(description, result)
    progress.update(task, advance=1)
    await q.put((name, taxid))


def update_otu(taxid: str, path: str):
    """
    Updates a otu.json's taxid key with either a taxon id or None
    depending on whether an id was able to be retrieved.

    Parameters:
        taxid (str): A taxon id for a given OTU if found, else None
        path (str): Path to a given OTU

    """
    with open(os.path.join(path, "otu.json"), 'r+') as f:
        otu = json.load(f)
        otu["taxid"] = int(taxid) if taxid else taxid
        f.seek(0)
        json.dump(otu, f, indent=4)


async def get_name_from_path(path: str, force_update: bool) -> Union[str, None]:
    """
    Given a path to an OTU, returns the name of the OTU. If a taxon id already exists
    for the OTU then it returns None.

    Parameters:
        path (str): Path to a src database directory
        force_update (bool): Flag to indicate whether an OTU should be searched for regardless
                             if they already have a taxon id

    Returns:
        The name of an OTU if a taxon id doesn't already exist in its JSON file, else None
        If the force_update flag is used the function will always return the OTU name

    """
    async with aiofiles.open(os.path.join(path, "otu.json"), "r") as f:
        otu = json.loads(await f.read())

        try:
            if otu["taxid"] is not None and not force_update:
                return None
            else:
                return otu["name"]
        except KeyError:
            return otu["name"]


def run(src_path: str, force_update: bool):
    """
    Creates and runs the event loop to asynchronously find taxon ids for all OTU in a src directory

    Parameters:
        src_path (str): Path to a database src directory
        force_update (bool): Flag to indicate whether all OTUs should be searched for regardless
                             if they already have a taxon id
    """
    asyncio.run(taxid(src_path, force_update))
