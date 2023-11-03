import json
from pathlib import Path
from typing import Optional
from structlog import get_logger, BoundLogger
from urllib.error import HTTPError

from virtool_cli.utils.reference import (
    get_otu_paths,
    get_isolate_paths,
    get_sequence_paths,
)
from virtool_cli.utils.ncbi import fetch_taxonomy_rank, fetch_upstream_record_taxids
from virtool_cli.utils.storage import parse_sequence
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


def search_by_otu_id(otu_id: str, catalog_path: Path) -> Optional[Path]:
    """
    Searches records for a matching id and returns the first matching path in the accession records

    :param otu_id: Unique OTU ID string
    :param catalog_path: Path to an accession catalog directory
    """
    matches = list(catalog_path.glob(f"*--{otu_id}.json"))
    if matches:
        return matches[0]

    return None


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


def get_otu_accessions(otu_path: Path) -> list:
    """
    Gets all accessions from an OTU directory and returns a list

    :param otu_path: Path to an OTU directory under a reference directory
    :return: A list of all accessions under an OTU
    """
    accessions = []

    for isolate_path in get_isolate_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            sequence = json.loads(sequence_path.read_text())
            accessions.append(sequence["accession"])

    return accessions


async def get_sequence_metadata(sequence_path: Path) -> dict:
    """
    Gets the accession length and segment name from a sequence file
    and returns it in a dict

    :param sequence_path: Path to a sequence file
    :return: A dict containing the sequence accession, sequence length and segment name if present
    """
    sequence = await parse_sequence(sequence_path)

    sequence_metadata = {
        "accession": sequence["accession"],
        "length": len(sequence["sequence"]),
    }

    segment = sequence.get("segment", None)
    if segment is not None:
        sequence_metadata["segment"] = segment

    return sequence_metadata


async def get_otu_accessions_metadata(otu_path: Path) -> dict:
    """
    Returns sequence metadata for all sequences present under an OTU

    :param otu_path: Path to an OTU directory under a reference directory
    :return: An accession-keyed dict containing all constituent sequence metadata
    """
    # get length and segment metadata from sequences
    all_metadata = {}

    for isolate_path in get_isolate_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            sequence_metadata = await get_sequence_metadata(sequence_path)

            accession = sequence_metadata["accession"]

            all_metadata[accession] = sequence_metadata

    return all_metadata
