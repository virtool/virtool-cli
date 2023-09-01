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
    Return all paths to accession listings relevant to the source reference.
    Uses the unique Virtool OTU ID to match between reference metadata and catalog listing.

    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    :return: A list of paths to relevant listings in the accession catalog
    """    
    included_listings = []
    for path in get_otu_paths(src_path):
        [ _, otu_id ] = (path.name).split('--')

        included_listings.append(
            search_by_otu_id(otu_id, catalog_path)
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

def search_by_otu_id(otu_id: str, catalog_path: Path):
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

def get_sequence_metadata(sequence_path: Path) -> list:
    """
    Gets accessions and sequence lengths from an OTU directory and returns a list

    :param otu_path: Path to an OTU directory
    """
    with open(sequence_path, "r") as f:
        sequence = json.load(f)

    sequence_metadata = {
        'accession': sequence['accession'],
        'length': len(sequence['sequence'])
    }
    if segment := sequence.get('segment', None):
        if segment is not None:
            sequence_metadata['segment'] = segment

    return sequence_metadata

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
    Checks each accession on the accession listing and requests metadata
    for each associated taxonomy ID.
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
    
def get_required_parts(schema: list):
    required_parts = []
    for part in schema:
        if part['required']:
            required_parts.append(part['name'])
    
    return required_parts

def measure_monopartite(sequence_metadata: dict):
    sequence_lengths = []
    print(sequence_metadata)
    for accession in sequence_metadata:
        metadata = sequence_metadata[accession]
        seq_length = metadata['length']
        sequence_lengths.append(seq_length)
    
    average_length = sum(sequence_lengths)/len(sequence_lengths)

    return int(average_length)

def measure_multipartite(sequence_metadata, part_list):
    part_total_dict = { element: [] for index, element in enumerate(part_list) }
            
    for accession in sequence_metadata:
        metadata = sequence_metadata[accession]
        seq_length = metadata['length']
        
        segment_name = metadata.get('segment', None)
        if segment_name is None:
            continue
        if segment_name in part_total_dict.keys():
            part_total_dict[segment_name].append(seq_length)
    
    length_dict = {}
    for segment in part_total_dict:
        if not part_total_dict[segment]:
            length_dict[segment] = 0
            continue
        average_length = sum(part_total_dict[segment]) / len(part_total_dict[segment])
        length_dict[segment] = int(average_length)

    # return int(average_length)
    return length_dict