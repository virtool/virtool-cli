from pathlib import Path
from Bio import Entrez
import asyncio
import structlog

from virtool_cli.utils.logging import DEBUG_LOGGER, DEFAULT_LOGGER
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL
from virtool_cli.catalog.listings import parse_listing, update_listing
from virtool_cli.catalog.helpers import find_taxid_from_accessions

base_logger = structlog.get_logger()


def run(catalog_path: Path, debugging: bool = False):
    """
    CLI entry point for accession.repair.run()

    :param catalog_path: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(catalog=str(catalog_path))

    logger.info("Repairing catalog...")

    asyncio.run(repair_catalog(catalog_path))


async def repair_catalog(catalog_path: Path):
    """
    Runs repair functions on the accession catalog

    :param catalog_path: Path to an accession catalog directory
    """
    await fill_missing_taxids(catalog_path)

    await rename_listings(catalog_path)

    return


async def fill_missing_taxids(catalog_path: Path):
    """
    Iterate through unmatched listings and run taxon ID extraction function.
    If a valid taxon ID is found, write to the listing.

    :param catalog_path: Path to an accession catalog directory
    """
    logger = base_logger

    for listing_path in catalog_path.glob("none--*.json"):
        logger = logger.bind(listing_path=str(listing_path.relative_to(catalog_path)))

        extracted_taxid = await find_taxid_from_accessions(listing_path, logger)

        if extracted_taxid:
            logger.debug(f"Found taxon ID {extracted_taxid}")

            listing = await parse_listing(listing_path)
            listing["taxid"] = extracted_taxid

            try:
                await update_listing(listing, listing_path)
                logger.info("Wrote new taxon ID to path")
            except Exception as e:
                logger.exception(e)
                continue


async def rename_listings(catalog_path: Path):
    """
    Renames listings where the taxon ID or OTU ID in the listing data
    no longer matches the listing's filename.

    :param catalog_path: Catalog listing data in dictionary form
    """
    logger = base_logger

    for listing_path in catalog_path.glob("*.json"):
        logger = logger.bind(listing_path=str(listing_path.relative_to(catalog_path)))

        [taxid, otu_id] = listing_path.stem.split("--")

        listing = await parse_listing(listing_path)

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
