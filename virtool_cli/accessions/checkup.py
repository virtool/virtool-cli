from pathlib import Path
import json
import structlog
from logging import INFO, DEBUG

from virtool_cli.accessions.helpers import get_catalog_paths, split_pathname

LISTING_KEYS = set(["_id", "accessions", "name", "schema", "schema", "taxid"])

b_logger = structlog.get_logger()


def run(catalog: Path, debugging: bool = False):
    """
    Check catalog for outstanding issues

    :param catalog: Path to an accession catalog directory
    :param debugging: Debugging flag
    """
    logger = b_logger.bind(command='checkup')

    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))
    
    if not catalog.exists():
        logger.critical('Catalog directory not found at path', path=str(catalog))
        return
    
    if unassigned := search_by_taxid('none', catalog):
        logger.warning(
            'Found listings without assigned NCBI taxon IDs.', 
            test='no_taxid',
            unassigned=unassigned)

    if duplicates := find_duplicates(catalog):
        logger.warning(
            'Found non-unique taxon IDs in catalog.', 
            test='duplicate_taxid',
            taxids=duplicates)
    
    if missing_uids := find_missing_uids(catalog):
        logger.warning(
            'Found listings with missing uids',
            test='null_uid',
            taxids=missing_uids
        )
    
    for listing_path in get_catalog_paths(catalog):
        with open(listing_path, "r") as f:
            listing = json.load(f)

        if missing_keys := check_keys(listing):
            logger.warning(
                'Entry is missing keys',
                missing_keys=missing_keys,
                path=str(listing_path.relative_to(catalog))
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

def find_duplicates(catalog: Path):
    """
    """
    duplicated_taxids = set()

    for listing_path in get_catalog_paths(catalog):
        [ taxid, otu_id ] = split_pathname(listing_path)

        logger = b_logger.bind(
            path=str(listing_path.relative_to(catalog.parent)), 
            otu_id=otu_id
        )

        if taxid != 'none':
            matches = search_by_taxid(taxid, catalog)

            if len(matches) > 1:
                logger.debug(f"Duplicate found!", taxid=taxid)
                duplicated_taxids.add(taxid)
    
    return duplicated_taxids

def find_missing_uids(catalog: Path):
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

if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-plant-viruses'
    src_path = project_path / 'src'
    # catalog_path = project_path / '.cache/catalog'
    catalog_path = Path(REPO_DIR) / 'ref-accession-catalog/catalog'
    # catalog_path = Path('/Users/sygao/Development/UVic/Virtool/TestSets/cotton_dupe/.cache/catalog')

    run(catalog=catalog_path, debugging=False)
