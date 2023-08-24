from pathlib import Path
import json
# from Bio import Entrez, SeqIO
import asyncio
import aiofiles
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ncbi import fetch_accession_uids, NCBI_REQUEST_INTERVAL
from virtool_cli.accessions.helpers import get_catalog_paths


def run(catalog: Path, debugging: bool = False):
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    base_logger.info('Repairing catalog...')
    
    asyncio.run(repair_catalog(catalog))

async def repair_catalog(catalog: Path):
    """
    """
    queue = asyncio.Queue()
    
    fetcher = asyncio.create_task(fetch_loop(catalog, queue))

    await asyncio.gather(*[fetcher], return_exceptions=True)

    asyncio.create_task(
        writer_loop(catalog, queue)
    )

    await queue.join()

    return

async def fetch_loop(catalog: Path, queue: asyncio.Queue):
    """
    """
    change_dict = {}

    for listing_path in get_catalog_paths(catalog):
        logger = base_logger.bind(listing_path=str(listing_path.relative_to(catalog)))

        with open(listing_path, "r") as f:
            listing = json.load(f)

        new_accessions = await fetch_missing_uids(listing)
        if new_accessions:
            change_dict[listing_path] = new_accessions
            logger.debug('Missing UIDs filled')

            await queue.put({ 'path': listing_path, 'listing': listing_path})
            await asyncio.sleep(NCBI_REQUEST_INTERVAL)

async def writer_loop(catalog_path: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the update function

    :param src_path: Path to a given reference directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    while True:
        packet = await queue.get()
        path = packet['path']
        listing = packet['listing']
        [ taxid, otu_id ] = path.name.split('--')
        logger = base_logger.bind(
            otu_id=otu_id, taxid=taxid, listing_path=str(path.relative_to(catalog_path))
        )

        try:
            update_listing(listing, path)
        except Exception as e:
            base_logger.exception('e')
        
        logger.info(f"Wrote updated listing to file")

        await asyncio.sleep(0.1)
        queue.task_done()

async def fetch_missing_uids(listing: dict):
    """
    """    
    change_flag = False

    for accession_list in ['included', 'excluded']:
        fetchlist = []
        for accession in listing['accessions'][accession_list]:
            if listing['accessions'][accession_list][accession] is None:
                fetchlist.append(accession)
    
        if not fetchlist:
            continue
        
        try:
            fetched_entries = await fetch_accession_uids(fetchlist)
        except:
            continue

        if not fetched_entries:
            continue
        
        for accession in fetched_entries:
            if fetched_entries[accession] is not None:
                listing['accessions'][accession_list][accession] = fetched_entries[accession]
                change_flag = True
            
    if change_flag:
        return listing
    else:
        return {}

async def update_listing(data, path):
    try:
        async with aiofiles.open(path, "w") as f: 
            await f.write(
                json.dumps(data, indent=2, sort_keys=True)
            )
    except Exception as e:
        return e

if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-plant-viruses'
    src_path = project_path / 'src'
    catalog_path = Path(REPO_DIR) / 'ref-accession-catalog/catalog'

    run(catalog_path, debug)
