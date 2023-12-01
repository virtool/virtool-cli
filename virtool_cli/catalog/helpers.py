from pathlib import Path
from structlog import get_logger, BoundLogger
from urllib.error import HTTPError

from virtool_cli.catalog.catalog import search_by_otu_id
from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.utils.ncbi import fetch_taxonomy_rank, fetch_upstream_record_taxids
from virtool_cli.catalog.listings import parse_listing


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
        [_, otu_id] = path.name.split("--")

        listing_path = search_by_otu_id(otu_id, catalog_path)
        if listing_path is not None:
            included_listings.append(listing_path)

    return included_listings


async def find_taxid_from_accessions(
    listing_path: Path, logger: BoundLogger = get_logger()
) -> str:
    """
    Checks each accession on the accession listing and requests metadata
    for each associated taxonomy ID.

    :param listing_path: Path to a listing in an accession catalog directory
    :param logger: Entry point for an existing BoundLogger
    :return: Taxon ID as string
    """
    listing = await parse_listing(listing_path)

    accessions = listing["accessions"]["included"]
    logger.debug("Searching unique IDs:", uids=accessions)

    try:
        records = await fetch_upstream_record_taxids(fetch_list=accessions)
    except HTTPError:
        logger.error("Failed to link taxon ID with nucleotide records.")
        return ""

    if not records:
        logger.warning("No taxon IDs found", taxids=records)
        return ""

    otu_taxids = []
    for taxid in records:
        rank = await fetch_taxonomy_rank(taxid)
        if rank == "species":
            otu_taxids.append(taxid)

    if not otu_taxids:
        logger.warning("No taxon IDs found", taxids=records)
        return ""

    if len(otu_taxids) > 1:
        logger.warning("Found multiple taxon IDs in this OTU", taxids=otu_taxids)
        return ""

    taxid = otu_taxids.pop()
    return taxid
