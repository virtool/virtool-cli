from pathlib import Path
import json
from Bio import Entrez, SeqIO
import asyncio
import aiofiles
from typing import Optional
from structlog import BoundLogger
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ncbi import fetch_taxonomy_rank, NCBI_REQUEST_INTERVAL
from virtool_cli.accessions.helpers import (
    get_catalog_paths, fix_listing_path, update_listing
)

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
    logger = base_logger.bind(catalog=str(catalog))

    # queue = asyncio.Queue()
    
    # fetcher = asyncio.create_task(fetch_loop(catalog, queue))

    # await asyncio.gather(*[fetcher], return_exceptions=True)

    # asyncio.create_task(
    #     writer_loop(catalog, queue)
    # )

    # await queue.join()

    await fetch_missing_taxids(catalog, logger)

    await rename_listings(catalog, logger)

    return

# async def fetch_loop(catalog: Path, queue: asyncio.Queue):
#     """
#     """
#     change_dict = {}

#     for listing_path in get_catalog_paths(catalog):
#         logger = base_logger.bind(listing_path=str(listing_path.relative_to(catalog)))

#         with open(listing_path, "r") as f:
#             listing = json.load(f)

#         new_accessions = await fetch_missing_uids(listing)
#         if new_accessions:
#             change_dict[listing_path] = new_accessions
#             logger.debug('Missing UIDs filled')

#             await queue.put({ 'path': listing_path, 'listing': listing_path})
#             await asyncio.sleep(NCBI_REQUEST_INTERVAL)

# async def writer_loop(catalog_path: Path, queue: asyncio.Queue) -> None:
#     """
#     Pulls packet dicts from the queue and calls the update function

#     :param src_path: Path to a given reference directory
#     :param queue: Queue of parsed OTU data awaiting processing
#     """
#     while True:
#         packet = await queue.get()
#         path = packet['path']
#         listing = packet['listing']
#         [ taxid, otu_id ] = path.name.split('--')
#         logger = base_logger.bind(
#             otu_id=otu_id, taxid=taxid, listing_path=str(path.relative_to(catalog_path))
#         )

#         try:
#             update_listing(listing, path)
#         except Exception as e:
#             base_logger.exception('e')
        
#         logger.info(f"Wrote updated listing to file")

#         await asyncio.sleep(0.1)
#         queue.task_done()

async def fetch_missing_taxids(catalog: Path, logger: BoundLogger):
    """
    """
    for listing_path in catalog.glob('none--*.json'):
        logger = logger.bind(listing_path=str(listing_path.relative_to(catalog)))
        extracted_taxid = await find_taxid_from_accession(listing_path, logger)
        if extracted_taxid is not None:
            logger.debug(f'Found taxon ID {extracted_taxid}')
            
            with open(listing_path, "r") as f:
                listing = json.load(f)
            listing['taxid'] = extracted_taxid
            
            try:
                await update_listing(listing, listing_path)
                logger.info('Wrote new taxon ID to path')
            except Exception as e:
                logger.exception(e)
                continue

async def rename_listings(catalog: Path, logger: BoundLogger):
    """
    Renames listings with improper 
    """
    # Rename misnamed listings
    for listing_path in catalog.glob('*.json'):
        logger = logger.bind(listing_path=str(listing_path.relative_to(catalog)))
        [ taxid, otu_id ] = listing_path.stem.split('--')

        with open(listing_path, "r") as f:
            listing = json.load(f)

        if taxid != str(listing['taxid']):
            logger.warning(
                f"Misnamed listing found: {taxid} != { listing['taxid'] }")
            fix_listing_path(listing_path, listing['taxid'], otu_id)
            taxid = listing['taxid']
        
        # if otu_id != listing['_id']:
        #     logger.warning(
        #         f"Misnamed listing found: {otu_id} != { listing['_id'] }")
        #     fix_listing_path(listing_path, taxid, listing['_id'])

async def find_taxid_from_accession(
    listing_path: Path, logger: BoundLogger = base_logger
):
    """
    """
    with open(listing_path, "r") as f:
        listing = json.load(f)

    logger = base_logger.bind(otu_id = listing['_id'])

    indexed_accessions = list(listing['accessions']['included'])
        
    try:
        records = await fetch_upstream_record_taxids(indexed_accessions)
    except Exception as e:
        logger.exception (e)
        return

    if not records:
        logger.warning('No taxon IDs found', taxids=records)
        return None
    
    otu_taxids = []
    for taxid in records:
        rank = await fetch_taxonomy_rank(taxid)
        if rank == 'species':
            otu_taxids.append(taxid)
    
    if not otu_taxids:
        logger.warning('No species-rank taxon IDs found', taxids=records)
        return None
    
    if len(otu_taxids) > 1:
        logger.warning('Found multiple taxon IDs in this OTU', taxids=otu_taxids)
        return None
    else:
        taxid = otu_taxids.pop()
        return taxid

async def fetch_upstream_record_taxids(
    fetch_list: list, 
) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank
    
    :param fetch_list: List of accession numbers to fetch from GenBank
    :param logger: Structured logger

    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=fetch_list, rettype="docsum")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise e
    
    await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    taxids = set()

    for r in record:
        taxids.add(int(r.get('TaxId')))
    
    return taxids

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
    
    project_path = Path(REPO_DIR) / 'ref-mini' #'ref-plant-viruses'
    src_path = project_path / 'src'
    # catalog_path = Path(REPO_DIR) / 'ref-accession-catalog/catalog'
    catalog_path = project_path / '.cache/catalog'

    run(catalog_path, debug)
