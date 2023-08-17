import json
import asyncio
import aiofiles
from pathlib import Path
from Bio import Entrez
import structlog
import logging

import virtool_cli.utils.logging
from virtool_cli.utils.ref import get_otu_paths, parse_otu
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL

logger = structlog.get_logger()


def run(src_path: Path, force_update: bool = False, debugging: bool = False):
    """
    Creates and runs the event loop to asynchronously find taxon ids for all OTU in a src directory

    :param src_path: Path to a given reference directory
    :param force_update: Flag to force update all of the OTU's taxids in a reference
    :return:
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    asyncio.run(taxid(src_path, force_update))

async def taxid(src_path: Path, force_update: bool = False):
    """
    Asynchronously finds taxon ids for all OTU in a given src directory and
    writes them to each respective otu.json file.

    :param src_path: Path to a given reference directory
    :param force_update: Flag to force update all of the OTU's taxids in a reference
    """
    queue = asyncio.Queue()

    otu_paths = get_otu_paths(src_path)
    
    if force_update:
        included_paths = otu_paths
    
    else:
        included_paths = filter_unidentified(src_path)
    
    logger.debug('Fetching taxon IDs in following paths...', included_paths=included_paths)
    
    if not included_paths:
        logger.info('No OTUs with missing taxon IDs to fetch', n_otus=len(included_paths))
        return
    
    logger.info(f'Requesting taxon IDs for {len(included_paths)}/{len(otu_paths)} OTUs', n_otus=len(included_paths))

    # Requests and retrieves matching taxon IDs from NCBI Taxonomy
    # and pushes results to queue
    fetcher = asyncio.create_task(
        fetcher_loop(included_paths, queue, force_update))
    
    # Pulls taxon id from queue and writes to otu.json
    asyncio.create_task(
        writer_loop(src_path, queue)
    )

    await asyncio.gather(fetcher, return_exceptions=True)

    await queue.join()

    # logger.info(
    #     f"Retrieved {counter} taxids for {len(included_paths)} OTUs", 
    #     n_updated=counter, n_otus=len(included_paths)
    # )

async def fetcher_loop(
    included_paths: list, 
    queue: asyncio.Queue, 
    force_update: bool = False
):
    """
    Loops through a list of specified OTU directories, 
    parses the data and pushes the name and taxon ID to the queue

    :param included_paths: List of paths to relevant OTUs
    :param queue: Queue of parsed OTU data awaiting processing
    """

    for path in included_paths:
        otu_data = parse_otu(path)

        if not force_update:
            taxid = otu_data.get('taxid', None)
            if taxid is not None:
                continue

        # Search for taxon id
        otu_name = otu_data.get('name')

        taxid = await fetch_by_taxid(name=otu_name)

        await log_results(name=otu_name, taxid=taxid)
        
        await queue.put({ 'path': path, 'name': otu_name, 'taxid': taxid})

        await asyncio.sleep(NCBI_REQUEST_INTERVAL)


async def writer_loop(src_path: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the update function

    :param src_path: Path to a given reference directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    while True:
        packet = await queue.get()
        taxid = packet['taxid']
        path = packet['path']
        [ _, otu_id ] = path.name.split('--')

        try:
            update_otu(taxid, path)
        except Exception as e:
            logger.exception('e')
        
        logger.info(f"Wrote taxon id '{taxid}' to file", 
            otu_id=otu_id, path=str(path.relative_to(src_path)))

        await asyncio.sleep(0.1)
        queue.task_done()

def filter_unidentified(src_path: Path) -> list:
    """
    Return OTU paths without assigned taxon ids

    :param src: Path to a reference directory
    :return: List of OTUs without assigned taxon ids
    """
    included_otus = []
    
    for path in get_otu_paths(src_path):
        otu_data = parse_otu(path)

        taxid = otu_data.get('taxid', None)
        
        if taxid is None:
            included_otus.append(path)
        
        else:
            continue
    
    return included_otus

async def fetch_by_taxid(name: str) -> int:
    """
    Searches the NCBI taxonomy database for a given OTU name and fetches and returns
    its taxon id. If a taxon id was not found the function should return None.

    :param name: Name of a given OTU
    :return: Taxonomy id for the given OTU
    """
    handle = Entrez.esearch(db="taxonomy", term=name)
    record = Entrez.read(handle)
    handle.close()

    try:
        taxid = int(record["IdList"][0])
    except IndexError:
        taxid = None

    return taxid

async def log_results(name: str, taxid: int=None):
    """
    Logs results of the taxid fetch from the NCBI taxonomy database

    :param name: Name of a given OTU
    :param taxid: Taxid for a given OTU if found, else None
    """
    otu_log = logger.bind(name=name)
    if taxid:
        otu_log.info('Taxon ID found', taxid=taxid)
    else:
        otu_log.debug('Taxon ID not found', taxid=None)


def update_otu(taxid: int, path: Path):
    """
    Updates a otu.json's taxid key with either a taxon id or None
    depending on whether an id was able to be retrieved.

    :param taxid: Taxid for a given OTU if found, else None
    :param path: Path to a given OTU in the reference directory
    """
    try:
        with open(path / "otu.json", "r+") as f:
            otu = json.load(f)
            otu["taxid"] = taxid if taxid else taxid
            f.seek(0)
            json.dump(otu, f, indent=4)

    except Exception as e:
        return e

async def get_name_from_path(path: Path, force_update: bool) -> str:
    """
    Given a path to an OTU, returns the name of the OTU. 
    If a taxon id already exists for the OTU, returns an empty string.

    :param path: Path to a given OTU
    :param force_update: Boolean flag used to determine if the function 
        should return the name of the OTU even if it already has a taxid
    :return: The OTU's name if it doesn't have a taxid assigned or force_update is True, 
        else an empty string
    """
    async with aiofiles.open(path / "otu.json", "r") as f:
        otu = json.loads(await f.read())

        try:
            if otu["taxid"] is not None and not force_update:
                return ""
            else:
                return otu["name"]
        except KeyError:
            return otu["name"]