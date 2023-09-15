import json
from pathlib import Path
from typing import Optional
from structlog import BoundLogger

from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.utils.ncbi import (
    fetch_taxonomy_rank, fetch_upstream_record_taxids
)


def get_catalog_paths(catalog: Path) -> list:
    """
    Return a list of paths to accession listings contained in an accession catalog.

    :param catalog_path: Path to an accession catalog directory
    :return: A list of paths representing the contents of the accession catalog.
    """
    return list(catalog.glob('*--*.json'))

def filter_catalog(src_path: Path, catalog_path: Path) -> list:
    """
    Return all paths to accession listings relevant to the source reference.
    Uses the unique Virtool OTU ID to match between reference metadata and catalog listing.

    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    :return: A list of paths to relevant listings in the accession catalog
    """    
    included_listings = []
    for path in get_otu_paths(src_path):
        [ _, otu_id ] = (path.name).split('--')
        
        listing_path = search_by_otu_id(otu_id, catalog_path)
        if listing_path is not None:
            included_listings.append(listing_path)
    
    return included_listings

def search_by_otu_id(otu_id: str, catalog_path: Path) -> Optional[Path]:
    """
    Searches records for a matching id and returns the first matching path in the accession records

    :param otu_id: Unique OTU ID string
    :param catalog_path: Path to an accession catalog directory
    """
    matches = list(catalog_path.glob(f'*--{otu_id}.json'))
    if matches:
        return matches[0]
    else:
        return None
    
async def find_taxid_from_accessions(
    listing_path: Path, logger: BoundLogger
) -> str:
    """
    Checks each accession on the accession listing and requests metadata
    for each associated taxonomy ID.

    :param listing_path: Path to a listing in an accession catalog directory
    :return: Taxon ID as string
    """
    with open(listing_path, "r") as f:
        listing = json.load(f)

    accessions = listing['accessions']['included']
    logger.debug('Searching uids:', uids=accessions)
    
    records = await fetch_upstream_record_taxids(fetch_list=accessions)
    if not records:
        logger.warning('No taxon IDs found', taxids=records)
        return None
    
    otu_taxids = []
    for taxid in records:
        rank = await fetch_taxonomy_rank(taxid)
        if rank == 'species':
            otu_taxids.append(taxid)
    
    if not otu_taxids:
        logger.warning('No taxon IDs found', taxids=records)
        return None
    
    if len(otu_taxids) > 1:
        logger.warning('Found multiple taxon IDs in this OTU', taxids=otu_taxids)
        return None
    else:
        taxid = otu_taxids.pop()
        return taxid