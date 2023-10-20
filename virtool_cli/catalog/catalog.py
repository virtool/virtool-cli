from pathlib import Path
from virtool_cli.catalog.listings import parse_listing


def get_catalog_paths(catalog_path: Path) -> list:
    """
    Return a list of paths to accession listings contained in an accession catalog.

    :param catalog_path: Path to an accession catalog directory
    :return: A list of paths representing the contents of the accession catalog.
    """
    return list(catalog_path.glob("*--*.json"))


def read_catalog(catalog_path: Path) -> dict:
    """
    Returns all catalog listings as a single dictionary keyed by Virtool ID

    :param catalog_path: Path to an accession catalog directory
    :return: A dictionary representing the accession catalog
    """
    catalog = {}

    for listing_path in catalog_path.glob("*--*.json"):
        [_, otu_id] = listing_path.name.split("--")

        catalog[otu_id] = parse_listing(listing_path)

    return catalog


def get_taxid_lookup(catalog_path: Path) -> dict:
    """
    Takes a catalog path and returns a dictionary of "taxid: otu_id" key-value pairs

    :param catalog_path: Path to an accession catalog directory
    :return: A lookup table matching the NCBI Taxonomy ID
    """
    lookup_table = {}

    for listing_path in catalog_path.glob("*--*.json"):
        [taxid, otu_id] = listing_path.name.split("--")
        lookup_table[taxid] = otu_id

    return lookup_table
