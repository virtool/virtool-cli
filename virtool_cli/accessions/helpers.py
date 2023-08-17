import json
from pathlib import Path
from virtool_cli.utils.ref import get_otu_paths, get_sequence_paths

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
