from pathlib import Path

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

def get_listings(catalog):
    return [listing for listing in catalog.glob('*.json')]

def split_pathname(path: Path):
    """
    Split a filename formatted as 'taxid--otu_id.json' into (taxid, otu_id)

    :param path: Path to a listing in an accession catalog
    """
    [ taxid, otu_id ] = (path.stem).split('--')

    return taxid, otu_id