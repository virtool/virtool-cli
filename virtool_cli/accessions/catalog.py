from pathlib import Path
import json
import asyncio
import structlog
from logging import INFO, DEBUG

from virtool_cli.utils.ref import (
    parse_otu,
    get_otu_paths,
    get_sequence_paths
)
# from virtool_cli.utils.ncbi import fetch_taxid, fetch_primary_ids
from virtool_cli.accessions.helpers import get_otu_accessions
from virtool_cli.accessions.initialize import initialize, generate_listing, write_listing

base_logger = structlog.get_logger()


def run(src_path: Path, catalog_path: Path, debugging: bool = False):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    :param debugging: Debugging flag
    """

    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))
    
    if not catalog_path.exists():
        base_logger.info('Initializing .cache/catalog directory...', catalog=str(catalog_path))
        catalog_path.mkdir()

    get_catalog(src_path, catalog_path)

def get_catalog(src: Path, catalog: Path):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    """
    logger = base_logger.bind(catalog=str(catalog))

    catalog_paths = [record for record in catalog.glob('*.json')]
    
    if catalog_paths:
        logger.info(
            'Catalog is already populated. Checking listings for accuracy...', catalogued=True)

        updated_list = check_catalog(src_path=src, catalog_path=catalog)
        
        if updated_list:
            logger.info(f'Updated {len(updated_list)} cache entries', 
                updated=updated_list)
        else:
            logger.info(f'Catalog is up to date')

    else:
        logger.info(
            'Empty catalog, initializing listings...', 
            catalog=str(catalog), catalogued=False)

        generate_catalog(src, catalog)
    
    # logger.info('Accessions written to cache', cache=str(cache))

def generate_catalog(src_path: Path, catalog_path: Path):
    """
    """
    asyncio.run(initialize(src_path, catalog_path))

def check_catalog(
    src_path: Path, catalog_path: Path
) -> list:
    """
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    changed_listings = []

    for otu_path in get_otu_paths(src_path):
        logger = base_logger.bind(
            otu_id=otu_path.name.split('--')[1],
            otu_path=str(otu_path.relative_to(otu_path.parents[1]))
        )

        otu_data = parse_otu(otu_path)
        taxid = otu_data.get('taxid', None)

        if taxid is None:
            logger.debug('taxid=null', name=otu_data['name'])
            continue

        listing_path = catalog_path / f"{taxid}--{otu_data.get('_id')}.json"
        
        logger = logger.bind(
            name=otu_data['name'], 
            taxid=taxid, 
            # listing_path=str(listing_path.relative_to(catalog_path.parent))
        )

        is_changed = check_listing(
            taxid=taxid,
            ref_accessions=get_otu_accessions(otu_path), 
            otu_path=otu_path, 
            listing_path=listing_path,
            logger=logger
        )
        if is_changed: changed_listings.append(taxid)
    
    return changed_listings

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