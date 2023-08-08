from pathlib import Path

def get_listings(catalog):
    return [listing for listing in catalog.glob('*.json')]

def split_pathname(path: Path):
    """
    Split a filename formatted as 'taxid--otu_id.json' into (taxid, otu_id)

    :param path: Path to a listing in an accession catalog
    """
    [ taxid, otu_id ] = (path.stem).split('--')

    return taxid, otu_id