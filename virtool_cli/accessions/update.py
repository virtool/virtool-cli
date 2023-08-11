from pathlib import Path
import json
import asyncio
import structlog
from logging import INFO, DEBUG

from virtool_cli.utils.ref import parse_otu, get_otu_paths, search_otu_by_id
from virtool_cli.utils.ncbi import fetch_accession_uids
from virtool_cli.accessions.helpers import parse_listing, split_pathname, get_otu_accessions, get_catalog_paths
from virtool_cli.accessions.initialize import generate_listing, write_listing

base_logger = structlog.get_logger()

def run(src: Path, catalog: Path):
    logger = base_logger.bind(task='update', catalog=str(catalog))
    
    logger.info(f"Starting catalog updater...")
    
    asyncio.run(update(src, catalog))

async def update(src: Path, catalog: Path):
    """
    """
    logger = base_logger.bind(task='update', catalog=str(catalog))

    check_catalog(src, catalog, logger)

    queue = asyncio.Queue()

    fetcher = asyncio.create_task(
        fetcher_loop(src, catalog, queue))
    
    asyncio.create_task(
        writer_loop(catalog, queue)
    )
    
    await asyncio.gather(fetcher, return_exceptions=True)

    await queue.join()

    logger.info(f"Catalog up to date")

async def fetcher_loop(src, catalog, queue):
    """
    """
    update_count = 0

    for listing_path in get_catalog_paths(catalog):
        [ taxid, otu_id ] = split_pathname(listing_path)

        logger = base_logger.bind(
            listing=str(listing_path.relative_to(listing_path.parents[1]))
        )

        otu_path = search_otu_by_id(src, otu_id)
        existing_accessions = set(get_otu_accessions(otu_path))

        listing_data = parse_listing(listing_path)
        listed_accessions = set(listing_data['accessions']['included'].keys())
        
        logger = logger.bind(
            otu_name=listing_data.get('name', ''),
            otu_id=otu_id,
            taxid=taxid
        )

        missing_accessions = existing_accessions.difference(listed_accessions)
        
        if listed_accessions != existing_accessions:
            logger.debug(f'difference = {missing_accessions}')

            indexed_accessions = listing_data['accessions']['included']

            for accession in missing_accessions:
                if accession in indexed_accessions.keys():
                    indexed_accessions.pop(accession)
                else:
                    indexed_accessions[accession] = None
                
                # logger.debug(f'new_index={indexed_accessions}')
                
                try:
                    missing_accessions = [key for (key, value) in indexed_accessions.items() if value is None]
                    # logger.debug(f'fetchlist={missing_accessions}')
                except Exception as e:
                    logger.exception(e)
                    return

                if missing_accessions:
                    new_indexed_accessions = await fetch_accession_uids(missing_accessions)
                    for accession in new_indexed_accessions:
                        indexed_accessions[accession] = new_indexed_accessions[accession]
                
                # listing_data['accessions']['included'] = indexed_accessions

            await queue.put({ 'path': listing_path, 'listing': listing_data } )
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
    while True:
        packet = await queue.get()

        with open(packet['path'], "w") as f:
            json.dump(packet['listing'], f, indent=2, sort_keys=True)

        await asyncio.sleep(0.1)
        queue.task_done()

def check_catalog(
    src_path: Path, catalog_path: Path,
    logger: structlog.BoundLogger = base_logger
) -> list:
    """
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    try:
        otu_ids = set([ (otu_path.name).split('--')[1] for otu_path in get_otu_paths(src)])
        catalog_ids = set([ split_pathname(listing_path)[1] for listing_path in get_catalog_paths(catalog)])
    except Exception as e:
        logger.exception(e)

    if otu_ids.issubset(catalog_ids):
        logger.info("All OTUs present in catalog")
    else:
        logger.warning("Missing OTUs", missing_otus=otu_ids.difference(catalog_ids))
        return

    # changed_listings = []

    # for otu_path in get_otu_paths(src_path):
    #     logger = base_logger.bind(
    #         otu_id=otu_path.name.split('--')[1],
    #         otu_path=str(otu_path.relative_to(otu_path.parents[1]))
    #     )

    #     otu_data = parse_otu(otu_path)
    #     taxid = otu_data.get('taxid', None)

    #     if taxid is None:
    #         logger.debug('taxid=null', name=otu_data['name'])
    #         continue

    #     listing_path = catalog_path / f"{taxid}--{otu_data.get('_id')}.json"
        
    #     logger = logger.bind(
    #         name=otu_data['name'], 
    #         taxid=taxid, 
    #         listing_path=str(listing_path.relative_to(catalog_path.parent))
    #     )

    #     is_changed = check_listing(
    #         taxid=taxid,
    #         ref_accessions=get_otu_accessions(otu_path), 
    #         otu_path=otu_path, 
    #         listing_path=listing_path,
    #         logger=logger
    #     )
    #     if is_changed: changed_listings.append(taxid)
    
    # return changed_listings

def check_listing(
    taxid: int, 
    ref_accessions: list, 
    otu_path: Path, 
    listing_path: Path,
    logger: structlog.BoundLogger
) -> bool:
    """
    :param taxid: Path to a reference directory
    :param ref_accessions: List of included accessions
    :param otu_path: Path to a OTU directory in a src directory
    :param listing_path: Path to a listing file in an accession catalog directory
    :param logger: Structured logger
    :return: True if changes have been made to the accession list, else False
    """
    if listing_path.exists():
        with open(listing_path, "r") as f:
            listing = json.load(f)


        ref_ids = set(ref_accessions)
        acc_ids = set(listing['accessions']['included'].keys())

        different_keys = ref_ids.difference(acc_ids)

        if len(different_keys) < 0:
            logger.debug('Changes found', different_keys=different_keys)

            listing.get('accessions')['included'] = ref_accessions
            
            update_listing(listing_path, ref_accessions, listing)
            logger.debug('Wrote new list')
            return True
        
        else:
            logger.debug('No changes found')

    else:
        logger.info(f'No accession record for {taxid}. Creating {listing_path.name}...')
        
        asyncio.run(
            initialize_listing(taxid, otu_path, ref_accessions, listing_path, logger)
        )
        
        return True
    
    return False

async def initialize_listing(
    taxid, otu_path, ref_accessions, listing_path, logger
):
    """
    """
    otu_data = parse_otu(otu_path)
    # otu_id = otu_data['_id']

    # accessions = get_otu_accessions(otu_path)

    new_record = await generate_listing(
        otu_data=otu_data, accession_list=ref_accessions)
    
    if not new_record:
        logger.error('Could not generate a listing for this OTU.')
        return
    
    await write_listing(taxid, new_record, catalog_path=listing_path.parent)

    return

def update_listing(
        path: Path,
        accessions: list,
        listing: dict, 
        indent: bool = True
):
    """
    Write accession file to cache directory

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