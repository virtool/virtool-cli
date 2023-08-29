from pathlib import Path
import json
from Bio import Entrez, SeqIO
import asyncio
import aiofiles
from structlog import BoundLogger
import logging
from urllib.error import HTTPError

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ncbi import (
    NCBI_REQUEST_INTERVAL,
    fetch_accession_uids, fetch_docsums, fetch_taxonomy_rank
)
from virtool_cli.accessions.helpers import (
    get_catalog_paths, fix_listing_path, update_listing, find_taxid_from_accession
)

def run(catalog: Path, debugging: bool = False):
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    base_logger.info('Fetching taxon IDs...')
    
    # wipe_taxids(catalog)
    asyncio.run(fetch_taxids(catalog))

async def fetch_taxids(catalog: Path):
    queue = asyncio.Queue()
    
    fetcher = asyncio.create_task(fetch_loop(catalog, queue))

    asyncio.create_task(writer_loop(catalog, queue))
    
    await asyncio.gather(*[fetcher], return_exceptions=True)

    await queue.join()

    return

async def fetch_loop(catalog: Path, queue: asyncio.Queue):
    """
    """
    base_logger.debug("Starting fetcher...", n_catalog=len(get_catalog_paths(catalog)))

    for listing_path in get_catalog_paths(catalog):
        logger = base_logger.bind(
            listing_path=str(listing_path.relative_to(catalog))
        )
        # logger.debug('Fetching taxon IDs from included accessions')
        try:
            extracted_taxid = await find_taxid_from_accession(listing_path, logger)
        except Exception as e:
            logger.exception(e)
            await asyncio.sleep(NCBI_REQUEST_INTERVAL)
            continue

        if extracted_taxid is not None:
            extracted_taxid = int(extracted_taxid)
            logger.debug(f'Found taxon ID: {extracted_taxid}')
            
            await queue.put({ 'path': listing_path, 'taxid': extracted_taxid })
        
        await asyncio.sleep(NCBI_REQUEST_INTERVAL)

async def writer_loop(catalog_path: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the update function

    :param src_path: Path to a given reference directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    write_logger = base_logger
    write_logger.debug("Starting writer...")
    while True:
        packet = await queue.get()
        path = packet['path']
        fetched_taxid = packet['taxid']
        [ old_taxid, otu_id ] = path.name.split('--')
        old_taxid = int(old_taxid)
        logger = base_logger.bind(
            otu_id=otu_id, listing_path=str(path.relative_to(catalog_path))
        )
        logger.debug(f'{packet}')

        if fetched_taxid != old_taxid:
            logger.error(f"All accessions are under taxon ID={fetched_taxid}, but current taxon ID is {old_taxid}")
        
        with open(path, "r") as f:
            listing = json.load(f)
        
        listing['taxid'] = fetched_taxid

        try:
            await update_listing(listing, path)
            logger.info(f"Wrote updated listing to file")
        except Exception as e:
            base_logger.exception(e)
        

        await asyncio.sleep(0.1)
        queue.task_done()

def wipe_taxids(catalog):
    """
    """
    for listing_path in catalog.glob('*--*.json'):

        with open(listing_path, "r") as f:
            listing = json.load(f)

        listing['taxid'] = "none"

        with open(listing_path, "w") as f: 
            f.write(json.dumps(listing, indent=2, sort_keys=True))

if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-plant-viruses'
    src_path = project_path / 'src'
    catalog_path = Path(REPO_DIR) / 'ref-accession-catalog/catalog'

    run(catalog_path, debug)