from pathlib import Path
import json
import asyncio
import structlog
from logging import INFO, DEBUG

from virtool_cli.utils.ref import parse_otu, get_otu_paths
from virtool_cli.utils.ncbi import fetch_taxid, fetch_accession_uids
from virtool_cli.accessions.helpers import get_otu_accessions

base_logger = structlog.get_logger()

def run(src: Path, catalog: Path):
    asyncio.run(initialize(src, catalog))

async def initialize(src: Path, catalog: Path):
    """
    """
    if not catalog.exists():
        catalog.mkdir()

    base_logger.debug(f"Starting catalog generation")

    queue = asyncio.Queue()

    fetcher = asyncio.create_task(
        fetcher_loop(src, queue))
    
    asyncio.create_task(
        writer_loop(catalog, queue)
    )
    
    await asyncio.gather(fetcher, return_exceptions=True)

    await queue.join()

    base_logger.info(f"Catalog generated")


async def fetcher_loop(src, queue):
    """
    """
    for otu_path in get_otu_paths(src):
        logger = base_logger.bind(
            otu_path=str(otu_path.relative_to(otu_path.parents[1]))
        )

        otu_data = parse_otu(otu_path)
        otu_id = otu_data['_id']
        
        logger = logger.bind(
            otu_name=otu_data.get('name', ''),
            otu_id=otu_id,
        )

        accessions = get_otu_accessions(otu_path)

        listing = await generate_listing(otu_data, accessions, logger)

        await queue.put({ 'otu_id': otu_id, 'listing': listing } )

async def writer_loop(catalog: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the update function

    :param src_path: Path to a given reference directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    while True:
        packet = await queue.get()

        otu_id = packet['otu_id']
        listing = packet['listing']
        taxid = listing['taxid']

        await write_listing(taxid, listing, catalog)

        await asyncio.sleep(0.1)
        queue.task_done()

async def generate_listing(
    otu_data: dict, 
    accession_list: list, 
    logger: structlog.BoundLogger = base_logger
) -> dict:
    """
    Generates a new listing for a given OTU and returns it as a dict

    :param otu_data: OTU data in dict form
    :param accession_list: list of included accesssions
    """
    catalog_listing = {}
    
    otu_id = otu_data.get('_id')
    taxid = otu_data.get('taxid', None)

    catalog_listing['_id'] = otu_id

    # Attempts to fetch the taxon id if none is found in the OTU metadata
    if taxid is None:
        logger.info('Taxon ID not found. Attempting to fetch from NCBI Taxonomy...')
        taxid = await fetch_taxid(otu_data.get('name', None))

        if taxid is None:
            catalog_listing['taxid'] = 'none'
            logger.info(f'Taxon ID not found. Setting taxid={taxid}')
            
        else:
            catalog_listing['taxid'] = int(taxid)
            logger.info(f'Taxon ID found. Setting taxid={taxid}')
    else:
        catalog_listing['taxid'] = int(taxid)

    catalog_listing['name'] = otu_data.get('name')
    
    schema = []
    for part in otu_data.get('schema'):
        if part.get('required'):
            schema.append(part.get('name'))
    catalog_listing['schema'] = schema
    
    # if len(schema) > 1:
    #     catalog_listing['multipartite'] = True
    # else:
    #     catalog_listing['multipartite'] = False
    
    try:
        indexed_accessions = await fetch_accession_uids(accession_list)
    except FileNotFoundError:
        return {}
    except RuntimeError as e:
        logger.exception(e)

    catalog_listing['accessions'] = {}
    catalog_listing['accessions']['included'] = indexed_accessions
    
    return catalog_listing

async def write_listing(
        taxid: int, 
        listing: dict, 
        catalog_path: Path,
        indent: bool = True):
    """
    Write accession file to cache directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param catalog: Path to an accession catalog
    :param indent: Indent flag
    """

    output_path = catalog_path / f"{taxid}--{listing['_id']}.json"
    logger = base_logger.bind(
        taxid=taxid, 
        listing_path=str(output_path.relative_to(output_path.parent)))

    logger.debug('Writing accession listing...')

    listing['accessions']['excluded'] = {}

    with open(output_path, "w") as f:
        json.dump(listing, f, indent=2 if indent else None, sort_keys=True)