import json
import aiofiles
from pathlib import Path
from typing import Optional
from structlog import BoundLogger

from virtool_cli.utils.ref import get_otu_paths, get_sequence_paths
from virtool_cli.utils.ncbi import fetch_taxonomy_rank, fetch_upstream_record_taxids

def get_catalog_paths(catalog) -> list:
    """
    Return a list of paths to accession listings contained in an accession catalog.

    :param catalog_path: Path to an accession catalog directory
    :return: A list of paths representing the contents of the accession catalog.
    """
    return list(catalog.glob('*.json'))

def filter_catalog(src_path, catalog_path) -> list:
    """
    Return paths for cached accession catalogues that are included in the source reference

    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    :return: A list of paths to relevant listings in the accession catalog
    """
    otu_paths = get_otu_paths(src_path)
    included_listings = []
    
    for path in otu_paths:
        [ _, otu_id ] = (path.name).split('--')

        included_listings.append(
            search_by_id(otu_id, catalog_path)
        )
    
    return included_listings

def parse_listing(path: Path):
    """
    Parse and return an accession listing.

    :param path: Path to a listing in an accession catalog
    """
    with open(path, "r") as f:
        listing = json.load(f)
    return listing

def split_pathname(path: Path):
    """
    Split a filename formatted as 'taxid--otu_id.json' into (taxid, otu_id)

    :param path: Path to a listing in an accession catalog
    """
    [ taxid, otu_id ] = (path.stem).split('--')

    return taxid, otu_id

def search_by_id(otu_id: str, catalog_path: Path):
    """
    Searches records for a matching id and returns the first matching path in the accession records

    :param otu_id: Unique OTU ID string
    :param catalog_path: Path to an accession catalog directory
    """
    matches = [listing for listing in catalog_path.glob(f'*--{otu_id}.json')]
    if matches:
        return matches[0]
    else:
        return FileNotFoundError
    
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

def fix_listing_path(path: Path, taxon_id: int, otu_id: str) -> Optional[str]:
    """
    Fixes each accession listing with the correct taxon ID as a label

    :param path: Path to a given reference directory
    :param otu: A deserialized otu.json
    """
    new_listing_name = f'{taxon_id}--{otu_id}.json'

    if path.name != new_listing_name:
        new_path = path.with_name(new_listing_name)
        path.rename(new_path)
        return new_path
    
async def update_listing(data, path):
    """
    """
    try:
        async with aiofiles.open(path, "w") as f: 
            await f.write(json.dumps(data, indent=2, sort_keys=True))
    except Exception as e:
        return e
    
async def find_taxid_from_accession(
    listing_path: Path, logger: BoundLogger
):
    """
    """
    with open(listing_path, "r") as f:
        listing = json.load(f)
    logger = logger.bind(otu_id = listing['_id'], listing=listing_path.name)

    indexed_uids = list(listing['accessions']['included'].values())
    logger.debug('Searching uids:', uids=indexed_uids)
    
    try:
        records = await fetch_upstream_record_taxids(indexed_uids)
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
        logger.warning('No taxon IDs found', taxids=records)
        return None
    
    if len(otu_taxids) > 1:
        logger.warning('Found multiple taxon IDs in this OTU', taxids=otu_taxids)
        return None
    else:
        taxid = otu_taxids.pop()
        return taxid