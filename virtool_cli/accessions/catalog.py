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
from virtool_cli.utils.ncbi import fetch_taxid, fetch_primary_ids

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

        all_current_accessions = generate_catalog(src)
    
        for otu_id in all_current_accessions:
            current_accessions = all_current_accessions.get(otu_id)

            logger.debug(
                f'Processing OTU _id {otu_id}', 
                taxid=current_accessions.get('taxid', None),
                current_accessions=current_accessions)
            
            try:
                write_listing(
                    taxid=current_accessions.get('taxid', None), 
                    record=current_accessions,
                    catalog_path=catalog
                )
            except Exception as e:
                logger.exception(e)
    
    
    # logger.info('Accessions written to cache', cache=str(cache))

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
        
        indexed_accessions = asyncio.run(fetch_primary_ids(ref_accessions))

        if set(indexed_accessions) == set(listing.get('accessions')['included']):
            logger.debug('No changes found')

        else:
            logger.debug('Changes found')
            # print(listing_path.name)
            # print(cached_record['accessions']['included'])
            # print(ref_accessions)

            listing.get('accessions')['included'] = ref_accessions
            
            update_listing(listing_path, ref_accessions, listing)
            logger.debug('Wrote new list')
            return True

    else:
        logger.info(f'No accession record for f{taxid}. Creating {listing_path.name}')
        
        new_record = generate_listing(
            otu_data=parse_otu(otu_path), 
            accession_list=ref_accessions)
        
        write_listing(taxid, new_record, catalog_path=listing_path.parent)
        
        return True
    
    return False

def generate_catalog(src: Path) -> dict:
    """
    Initialize an accession catalog by pulling accessions from all OTU directories

    :param src: Path to a reference directory
    :return: A dictionary of all accessions in each OTU, indexed by internal OTU id
    """
    catalog = {}

    for otu_path in get_otu_paths(src):
        logger = base_logger.bind(
            path=str(otu_path.relative_to(otu_path.parents[1]))
        )

        otu_data = parse_otu(otu_path)
        otu_id = otu_data['_id']
        
        logger = logger.bind(
            name=otu_data.get('name', ''),
            otu_id=otu_id,
        )
        
        accessions = get_otu_accessions(otu_path)

        catalog[otu_id] = generate_listing(otu_data, accessions, logger)
    
    return catalog

def generate_listing(
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
        taxid = asyncio.run(fetch_taxid(otu_data.get('name', None)))

        if taxid is None:
            catalog_listing['taxid'] = 'none'
            logger.info(f'Taxon ID not found. Setting taxid={taxid}')
            
        else:
            catalog_listing['taxid'] = int(taxid)
            logger.info(f'Taxon ID found. Setting taxid={taxid}')
    else:
        catalog_listing['taxid'] = int(taxid)

    catalog_listing['name'] = otu_data.get('name')

    schema = otu_data.get('schema', [])
    if len(schema) > 1:
        catalog_listing['multipartite'] = True
    else:
        catalog_listing['multipartite'] = False
    
    indexed_accessions = asyncio.run(fetch_primary_ids(accession_list))
    # print(indexed_accessions)

    catalog_listing['accessions'] = {}
    catalog_listing['accessions']['included'] = indexed_accessions
    
    return catalog_listing

def get_otu_accessions(otu_path: Path) -> list:
    """
    Gets all accessions from an OTU directory and returns a list

    :param otu_path: Path to an OTU directory
    """
    accessions = []
    
    for isolate_path in get_otu_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            with open(sequence_path, "r") as f:
                sequence = json.load(f)
            accessions.append(sequence['accession'])

    return accessions

def write_listing(
        taxid: int, 
        record: dict, 
        catalog_path: Path,
        indent: bool = True):
    """
    Write accession file to cache directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param catalog: Path to an accession catalog
    :param indent: Indent flag
    """

    logger = base_logger.bind(taxid=taxid)

    output_path = catalog_path / f"{taxid}--{record['_id']}.json"

    logger.debug('Writing accession', accession_path=output_path)

    record['accessions']['excluded'] = {}

    with open(output_path, "w") as f:
        json.dump(record, f, indent=2 if indent else None)
    
    return

def update_listing(
        path: Path,
        accessions: list,
        listing: dict, 
        indent: bool = True):
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