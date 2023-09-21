from pathlib import Path
import json
from Bio import Entrez
import asyncio
from structlog import BoundLogger
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL
from virtool_cli.accessions.listings import update_listing
from virtool_cli.accessions.helpers import find_taxid_from_accessions


def run(catalog_path: Path, debugging: bool = False):
    """
    CLI entry point for accession.repair.run()

    :param catalog_path: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    base_logger.info("Repairing catalog...")

    asyncio.run(repair_catalog(catalog_path))


async def repair_catalog(catalog_path: Path):
    """
    Runs repair functions on the accession catalog

    :param catalog_path: Path to an accession catalog directory
    """
    logger = base_logger.bind(catalog=str(catalog_path))

    await fill_missing_taxids(catalog_path, logger)

    await rename_listings(catalog_path, logger)

    return


async def fill_missing_taxids(catalog_path: Path, logger: BoundLogger = base_logger):
    """
    Iterate through unmatched listings and run taxon ID extraction function.
    If a valid taxon ID is found, write to the listing.

    :param listing: Catalog listing data in dictionary form
    :param logger: Optional entry point for a shared BoundLogger
    """
    for listing_path in catalog_path.glob("none--*.json"):
        logger = logger.bind(listing_path=str(listing_path.relative_to(catalog_path)))

        extracted_taxid = await find_taxid_from_accessions(listing_path, logger)

        if extracted_taxid is not None:
            logger.debug(f"Found taxon ID {extracted_taxid}")

            with open(listing_path, "r") as f:
                listing = json.load(f)
            listing["taxid"] = extracted_taxid

            try:
                await update_listing(listing, listing_path)
                logger.info("Wrote new taxon ID to path")
            except Exception as e:
                logger.exception(e)
                continue


async def rename_listings(catalog_path: Path, logger: BoundLogger = base_logger):
    """
    Renames listings where the taxon ID or OTU ID in the listing data
    no longer matches the listing's filename.

    :param catalog_path: Catalog listing data in dictionary form
    :param logger: Optional entry point for a shared BoundLogger
    """
    for listing_path in catalog_path.glob("*.json"):
        logger = logger.bind(listing_path=str(listing_path.relative_to(catalog_path)))

        [taxid, otu_id] = listing_path.stem.split("--")

        listing = json.loads(listing_path.read_text())

        if taxid != str(listing["taxid"]) or otu_id != listing["_id"]:
            old_name = listing_path.name

            if listing["taxid"] < 0:
                new_path = listing_path.with_name(f"none--{otu_id}.json")
            else:
                new_path = listing_path.with_name(f"{taxid}--{otu_id}.json")

            listing_path.rename(new_path)

            logger.info(f"Renamed {old_name} to {new_path.name}", listing=new_path.name)


async def fetch_upstream_record_taxids(fetch_list: list) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank

    :param fetch_list: List of accession numbers to fetch from GenBank
    :return: A list of GenBank data converted from XML to dicts if possible,
        else an empty list
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=fetch_list, rettype="docsum")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise e

    await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    taxids = []

    for r in record:
        taxids.append(int(r.get("TaxId")))

    return taxids
