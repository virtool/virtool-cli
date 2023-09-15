from pathlib import Path
import json
import asyncio
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL
from virtool_cli.accessions.helpers import (
    get_catalog_paths, update_listing, find_taxid_from_accessions
)


def run(catalog: Path, debugging: bool = False):
    """
    CLI entry point for accession.taxid.run()

    :param catalog_path: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    base_logger.info('Fetching taxon IDs...')
    
    # wipe_taxids(catalog)
    asyncio.run(fetch_taxids(catalog))

async def fetch_taxids(catalog: Path):
    """
    Extracts taxon IDs from the included accessions.
    Runs the taxon ID request fetcher and the listing writer. 

    :param catalog_path: Path to an accession catalog directory
    """
    queue = asyncio.Queue()
    
    fetcher = asyncio.create_task(fetcher_loop(catalog, queue))

    asyncio.create_task(writer_loop(catalog, queue))
    
    await asyncio.gather(*[fetcher], return_exceptions=True)

    await queue.join()

async def fetcher_loop(catalog: Path, queue: asyncio.Queue):
    """
    Iterates through all listings in a catalog directory, reads the included accessions 
    and requests the taxon ID data from NCBI

    :param catalog_path: Path to a catalog directory
    :param queue: Queue containing the listing path and all taxon IDs extracted from the included accessions
    """
    base_logger.debug("Starting fetcher...", n_catalog=len(get_catalog_paths(catalog)))

    for listing_path in catalog.glob('*.json'):
        logger = base_logger.bind(
            listing_path=str(listing_path.relative_to(catalog))
        )

        try:
            extracted_taxid = await find_taxid_from_accessions(
                listing_path, logger)
            
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
    Pulls packet dicts from the queue and updates the listing file accordingly

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
            logger.error(
                f"All accessions are under taxon ID={fetched_taxid}," + \
                f"but current taxon ID is {old_taxid}"
            )
        
        with open(path, "r") as f:
            listing = json.load(f)
        
        listing['taxid'] = fetched_taxid

        try:
            await update_listing(listing, path)
            logger.info("Wrote updated listing to file")
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