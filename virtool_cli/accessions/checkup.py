from pathlib import Path
import json
import structlog
import logging
# from logging import INFO, DEBUG

from virtool_cli.utils.logging import base_logger
from virtool_cli.accessions.helpers import get_catalog_paths, split_pathname

LISTING_KEYS = set(["_id", "accessions", "name", "schema", "schema", "taxid"])

def run(catalog: Path, debugging: bool = False):
    """
    Check catalog for outstanding issues. 
    Can also be used to diagnose issues with the corresponding reference.

    :param catalog: Path to an accession catalog directory
    :param debugging: Debugging flag
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    logger = base_logger.bind(catalog=str(catalog))
    
    if not catalog.exists():
        logger.critical('Catalog directory not found at path', path=str(catalog))
        return
    
    if unassigned := search_by_taxid('none', catalog):
        logger.debug(
            'Found listings without assigned NCBI taxon IDs.', 
            test='no_taxid',
            unassigned=unassigned)

    if duplicate_otus := find_duplicate_otus(catalog, logger):
        logger.warning(
            'Found non-unique taxon IDs in catalog.', 
            test='duplicate_taxid',
            taxids=duplicate_otus)
    
    if missing_uids := find_missing_uids(catalog, logger):
        logger.error(
            'Found listings with missing or malformed uids',
            test='null_uid',
            taxids=missing_uids
        )

    if duplicate_accessions := find_duplicate_accessions(catalog, logger):
        logger.error(
            'Found non-unique accessions in listings', 
            test='duplicate_accessions',
            otus=duplicate_accessions)
    
    for listing_path in get_catalog_paths(catalog):
        with open(listing_path, "r") as f:
            listing = json.load(f)

        if missing_keys := check_keys(listing):
            logger.warning(
                'Entry is missing keys',
                missing_keys=missing_keys,
                path=str(listing_path.relative_to(catalog, logger))
            )
        
        if not check_schema(listing):
            logger.warning('Schema is empty')

def check_keys(listing: dict):
    """
    """
    if not LISTING_KEYS.issubset(listing.keys()):
        missing_keys = LISTING_KEYS.difference(listing.keys())
        return list(missing_keys)
    else:
        return []

def check_schema(listing: dict):
    """
    """
    if type(listing['schema']) != list:
        return False
    
    if not listing:
        return False

    return True


def search_by_taxid(taxid, catalog_path: Path) -> list:
    """
    Searches records for a matching taxon id and returns all matching paths in the accession records

    :param otu_id: Unique taxon ID
    :param catalog_path: Path to an accession catalog directory
    """
    
    matches = [str(listing.relative_to(catalog_path)) 
        for listing in catalog_path.glob(f'{taxid}--*.json')]
    
    if matches:
        return matches
    else:
        return []

def find_duplicate_otus(catalog: Path, logger = base_logger):
    """
    """
    duplicated_taxids = set()

    for listing_path in get_catalog_paths(catalog):
        [ taxid, otu_id ] = split_pathname(listing_path)

        logger = logger.bind(
            path=str(listing_path.relative_to(catalog.parent)), 
            otu_id=otu_id
        )

        if taxid != 'none':
            matches = search_by_taxid(taxid, catalog)

            if len(matches) > 1:
                logger.debug(f"Duplicate taxon id found!", taxid=taxid)
                duplicated_taxids.add(taxid)
    
    return duplicated_taxids

def find_missing_uids(catalog: Path, logger = base_logger):
    """
    """
    missing_uids = set()
    for listing_path in get_catalog_paths(catalog):
        relative_path = str(listing_path.relative_to(listing_path.parent))
        with open(listing_path, "r") as f:
            listing = json.load(f)
        
        if type(listing['accessions']['included']) is not dict:
            missing_uids.add(relative_path)

        if type(listing['accessions']['excluded']) is not dict:
            missing_uids.add(relative_path)

        for accession in listing['accessions']['included']:
            if listing['accessions']['included'][accession] is None:
                missing_uids.add(relative_path)

        for accession in listing['accessions']['excluded']:
            if listing['accessions']['excluded'][accession] is None:
                missing_uids.add(relative_path)

            
    
    return missing_uids

def find_duplicate_accessions(catalog: Path, logger = base_logger):
    """
    """
    duplicate_accessions = []
    for listing_path in get_catalog_paths(catalog):
        with open(listing_path, "r") as f:
            listing = json.load(f)
        otu_id = listing['_id']
        logger = logger.bind(otu_id=listing['_id'], taxid=listing['taxid'], name=listing['name'])

        for alist_type in ['included', 'excluded']:
            keys = listing['accessions'][alist_type].keys()

            for key in keys:
                key_split = key.split('.')
                if len(key_split) > 1:
                    [ accession, version ] = key_split
                    if accession in keys:
                        logger.debug('Duplicate accession number found')
                        duplicate_accessions.append(otu_id)
    
    return duplicate_accessions
