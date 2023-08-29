from pathlib import Path
import json
import asyncio
import structlog
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import parse_otu, get_otu_paths, search_otu_by_id
from virtool_cli.accessions.initialize import generate_listing, write_listing
from virtool_cli.accessions.helpers import (
    parse_listing, split_pathname, get_otu_accessions, 
    get_catalog_paths, filter_catalog,
    # fix_listing_path
)

def run(src: Path, catalog: Path, debugging: bool = False):
    """
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    logger = base_logger.bind(task='update', catalog=str(catalog))
    
    logger.info(f"Starting catalog updater...")
    
    asyncio.run(update(src, catalog))

async def update(src: Path, catalog: Path):
    """
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    logger = base_logger.bind(task='update', catalog=str(catalog))

    await check_catalog(src, catalog, logger)

    queue = asyncio.Queue()

    fetcher = asyncio.create_task(
        fetcher_loop(src, catalog, queue))
    
    asyncio.create_task(
        writer_loop(catalog, queue))
    
    await asyncio.gather(fetcher, return_exceptions=True)

    await queue.join()

    logger.info(f"Catalog up to date")

async def fetcher_loop(src: Path, catalog: Path, queue: asyncio.Queue):
    """
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    fetch_logger = base_logger.bind(
        loop='fetcher',
        # catalog=str(catalog)
    )
    fetch_logger.debug("Starting fetcher...", catalog_size=len(get_catalog_paths(catalog)))

    update_count = 0

    for listing_path in filter_catalog(src, catalog):
        # fetch_logger.debug(f'{listing_path}')
        
        [ taxid, otu_id ] = split_pathname(listing_path)

        logger = fetch_logger.bind(listing=listing_path.name)
        
        otu_path = search_otu_by_id(src, otu_id)
        existing_accessions = set(get_otu_accessions(otu_path))

        listing = parse_listing(listing_path)
        logger = logger.bind(
            otu_name=listing['name'],
            otu_id=otu_id,
            taxid=taxid
        )

        listed_accessions = set(listing['accessions']['included'])
        
        not_in_reference = listed_accessions.difference(existing_accessions)
        not_in_listing = existing_accessions.difference(listed_accessions)
        indexed_accessions = listing['accessions']['included']
        # indexed_excluded = listing['accessions']['excluded']

        if not_in_reference:
            logger.debug(f'reference is missing accessions: {not_in_reference}')
            for accession in not_in_reference:
                indexed_accessions.pop(accession)

        if not_in_listing:
            logger.debug(f'listing is missing accessions: {not_in_listing}')

            for accession in not_in_listing:
                if accession in indexed_accessions:
                    logger.debug(f'Removing missing accession {accession}')
                    indexed_accessions.pop(accession)
                else:
                    indexed_accessions.append(accession)

                if not_in_listing:
                    for accession in not_in_listing:
                        indexed_accessions.append(accession)
                
                listing['accessions']['included'] = indexed_accessions

            await queue.put({ 'path': listing_path, 'listing': listing } )
            logger.debug(f'Pushed updated listing to queue', updated_accessions=indexed_accessions)
            update_count += 1

    base_logger.debug(
        f'Pushed {update_count} updates to writer', 
        n_updated=update_count, catalog_path=str(catalog))

async def writer_loop(catalog: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the update function

    :param src_path: Path to a given reference directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    write_logger = base_logger.bind(
        loop='writer',
        catalog=str(catalog)
    )
    write_logger.debug("Starting writer...")

    while True:
        packet = await queue.get()
        listing_path = packet['path']
        listing = packet['listing']
        write_logger.bind(
            listing_path=str(listing_path.relative_to(catalog))
        )

        write_logger.debug(f'New listing:\n{listing}')

        with open(packet['path'], "w") as f:
            json.dump(packet['listing'], f, indent=2, sort_keys=True)

        await asyncio.sleep(0.1)
        queue.task_done()

async def check_catalog(
    src_path: Path, catalog_path: Path,
    logger: structlog.BoundLogger = base_logger
) -> list:
    """
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    try:
        otu_ids = set([ (otu_path.name).split('--')[1] for otu_path in get_otu_paths(src_path)])
        catalog_ids = set([ split_pathname(listing_path)[1] for listing_path in get_catalog_paths(catalog_path)])
    except Exception as e:
        logger.exception(e)

    if otu_ids.issubset(catalog_ids):
        logger.info("All OTUs present in catalog")
        return True

    # If if we're missing OTUs, keep going
    missing_otus = otu_ids.difference(catalog_ids)
    logger.warning("Missing listings for OTUs", missing_otus=missing_otus)

    for otu_id in missing_otus:
        otu_path = [ path for path in src_path.glob(f'*--{otu_id}') ][0]
        await add_listing(otu_path, catalog_path, logger=logger)
    
    return False

async def add_listing(
    otu_path: Path, catalog: Path, 
    logger: structlog.BoundLogger = base_logger
):
    """
    If no listing is found for a particular OTU, generate and write a new listing
    to the catalog
    """
    otu_data = parse_otu(otu_path)
    
    logger.info(f"No accession record for {otu_data['taxid']}.")

    ref_accessions = get_otu_accessions(otu_path)

    new_record = await generate_listing(
        otu_data=otu_data, accession_list=ref_accessions)
    
    if not new_record:
        logger.error('Could not generate a listing for this OTU.')
        return
    
    await write_listing(otu_data['taxid'], new_record, catalog_path=catalog)

    return

def update_listing(
        path: Path,
        accessions: list,
        listing: dict, 
        indent: bool = True
):
    """
    Write accession file to catalog directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param listing: 
    :param indent: Indent flag
    """
    taxid_log = base_logger.bind(taxid=path.stem)

    taxid_log.debug('Updating accession', accession_path=path)

    listing['accessions']['included'] = accessions

    with open(path, "w") as f:
        json.dump(listing, f, indent=2 if indent else None)
    
    return