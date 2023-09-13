from pathlib import Path
import json
import asyncio
from structlog import BoundLogger
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ncbi import get_spelling
from virtool_cli.accessions.helpers import get_catalog_paths, search_by_taxid

LISTING_KEYS = set(["_id", "accessions", "name", "schema", "taxid"])


def run(catalog: Path, debugging: bool = False):
    """
    Entry point for CLI. Sets logging levels.

    :param catalog: Path to an accession catalog directory
    :param debugging: Debugging flag
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    run_tests(catalog)

def run_tests(catalog: Path):
    """
    Check catalog for outstanding issues. 
    Can also be used to diagnose issues with the corresponding reference.

    :param catalog: Path to an accession catalog directory
    """
    logger = base_logger.bind(catalog=str(catalog))
    
    logger.info('Checking for missing data within all listings...')
    check_missing_data(catalog, logger)

    logger.info('Checking for listings without taxon IDs...', test='no_taxid')
    if unassigned := search_by_taxid('none', catalog):
        logger.warning(
            'Found listings without assigned NCBI taxon IDs.', 
            test='no_taxid',
            unassigned=unassigned)

    logger.info('Checking for duplicate OTUs...', test='duplicate_taxid')
    if duplicate_otus := find_shared_taxids(catalog, logger):
        logger.warning(
            'Found non-unique taxon IDs in catalog.', 
            test='duplicate_taxid',
            taxids=duplicate_otus)

    logger.info(
        'Checking for listings containing duplicate accessions...', 
        test='duplicate_accessions')
    if duplicate_accessions := find_duplicate_accessions(catalog, logger):
        logger.warning(
            'Found non-unique accessions in listings', 
            test='duplicate_accessions',
            otus=duplicate_accessions)

    logger.info(
        'Requesting spelling suggestions for unmatchable listings...',
        test='suggest_spellings')
    asyncio.run(suggest_spellings(catalog))

def check_missing_data(
    catalog: Path, 
    logger: BoundLogger = base_logger
):
    """
    Checks all catalog listings for missing keys and schema data.

    :param catalog: Path to a catalog directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    for listing_path in get_catalog_paths(catalog):
        logger = logger.bind(listing=listing_path.stem)

        listing = json.loads(listing_path.read_text())

        if missing_keys := check_keys(listing):
            logger.warning(
                'Entry is missing keys',
                missing_keys=missing_keys,
                path=str(listing_path.relative_to(catalog, logger))
            )
        
        if not check_schema(listing):
            logger.warning('Schema field is empty.')

def check_keys(listing: dict):
    """
    Checks a listing to ensure that all necessary keys are present.
    Returns an empty list if all keys are present. 

    :param listing: Catalog listing data in dictionary form
    :returns: A list of all missing keys
    """
    if not LISTING_KEYS.issubset(listing.keys()):
        missing_keys = LISTING_KEYS.difference(listing.keys())
        return list(missing_keys)
    else:
        return []

def check_schema(listing: dict):
    """
    Checks a listing to ensure that the schema field has data in it and returns a boolean

    :param listing: Catalog listing data in dictionary form
    :returns: A boolean confirming the existence of a filled-out schema
    """
    if type(listing['schema']) != list:
        return False
    
    if not listing:
        return False

    return True


def find_shared_taxids(
    catalog_path: Path, 
    logger: BoundLogger = base_logger
) -> list:
    """
    Find OTUs that share a taxon ID and return them in a list

    :param catalog: Path to a catalog directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: A list of taxon IDs that have multiple OTUs under them
    """
    duplicated_taxids = set()

    for listing_path in catalog_path.glob('*--*.json'):
        [ taxid, otu_id ] = listing_path.stem('--')

        logger = logger.bind(
            path=str(listing_path.relative_to(catalog_path.parent)), 
            otu_id=otu_id
        )

        if taxid != 'none':
            matches = search_by_taxid(taxid, catalog_path)

            if len(matches) > 1:
                logger.debug("Duplicate taxon id found", taxid=taxid)
                duplicated_taxids.add(taxid)
    
    return list(duplicated_taxids)

def find_duplicate_accessions(catalog: Path, logger = base_logger):
    """
    Checks catalog listings for accessions that have been listed twice.
    Can be used to identify redundant sequences in the reference directory.

    :param catalog: Path to a catalog directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    duplicate_accessions = []

    for listing_path in catalog.glob('*--*.json'):
        listing = json.loads(listing_path.read_text())
            
        logger = logger.bind(
            listing=listing_path.name, name=listing['name'])
               
        logger.debug('Checking for non-unique accessions within included/excluded lists...')
        for alist_type in ['included', 'excluded']:
            accession_list = listing['accessions'][alist_type]
            accession_set = set(accession_list)

            if len(accession_set) < len(accession_list):
                logger.warning('Contains non-unique accessions')
                duplicate_accessions.append(listing_path.name)
        
        logger.debug('Checking included list against excluded list for eliminations...')
        for versioned_accession in listing['accessions']['included']:
            [ accession, version ] = versioned_accession.split('.')
            if accession in accession_set:
                logger.warning(
                    f"Included accession '{versioned_accession}' is on the exclusion list."
                )
                duplicate_accessions.append(listing_path.name)
    
    return duplicate_accessions

async def suggest_spellings(catalog: Path, logger = base_logger):
    """
    Evaluates the names of OTUs without retrievable taxon IDs
    and queries Entrez ESpell for alternatives.

    :param catalog: Path to a catalog directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    for listing_path in catalog.glob('none--*.json'):
        with open(listing_path, "r") as f:
            listing = json.load(f)
        
        logger = logger.bind(
            path = str(listing_path.relative_to(catalog)),
            otu_id=listing['_id'], 
            taxid=listing['taxid'], 
            current_name=listing['name'])
        
        current_name = listing['name']
        
        if '-' in current_name:
            current_name = current_name.split('-')[0]
        
        new_spelling = await get_spelling(current_name)

        if new_spelling != current_name.lower():
            logger.info(f'Try: {new_spelling}', potential_name=new_spelling)