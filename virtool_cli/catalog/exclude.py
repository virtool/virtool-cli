from pathlib import Path
import json
import asyncio
import aiofiles
from Bio import Entrez, SeqIO
import structlog
from structlog import BoundLogger
from urllib.error import HTTPError

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.catalog.listings import parse_listing

base_logger = structlog.get_logger()


def run(catalog_path: Path, debugging: bool = False):
    """
    CLI entry point for catalog.exclude.run()

    :param catalog_path: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(catalog=str(catalog_path), verbose=debugging)

    logger.info("Starting RefSeq filter...", catalog=str(catalog_path))

    asyncio.run(filter_refseq_accessions(catalog_path))


async def filter_refseq_accessions(catalog_path: Path):
    """
    Automatically excludes the former accession numbers of RefSeq sequences
    Runs RefSeq filter and exclusion on all listings in an accession catalog

    :param catalog_path: Path to a catalog directory
    """
    logger = base_logger
    for listing_path in catalog_path.glob("*--*.json"):
        logger = logger.bind(listing_path=listing_path.name)

        await filter_refseq_otu(listing_path, logger)


async def filter_refseq_otu(listing_path: Path, logger: BoundLogger = base_logger):
    """
    Parses a catalog listing, checks the included accession list
    for RefSeq-formatted accessions and requests GenBank for data.
    Checks the GenBank Comment data for a duplicate accession,
    then adds that accession to the excluded accession list.

    :param listing_path: Path to an individual listing in an accession catalog directory
    :param logger: Optional entry point for an existing BoundLogger
    """
    listing = await parse_listing(listing_path)

    logger = logger.bind(name=listing["name"])

    refseq_accessions = []
    for accession in listing["accessions"]["included"]:
        # Only RefSeq accessions include underscores
        if "_" in accession:
            refseq_accessions.append(accession)

    if not refseq_accessions:
        return

    logger.debug("RefSeq accessions found", refseq=refseq_accessions)
    try:
        included_records = await fetch_upstream_records(refseq_accessions)
    except Exception as e:
        logger.exception(e)
        return

    excluded_set = set(listing["accessions"]["excluded"])

    new_exclusions = []

    for record in included_records:
        accession = record["id"]
        redundant_accession = in_refseq(
            comments=record.annotations["comment"], excluded=excluded_set
        )

        if redundant_accession:
            logger.debug(
                f"Equivalent accession '{redundant_accession}' found for {accession}"
            )
            new_exclusions.append(redundant_accession)

    if new_exclusions:
        logger.info("Adding accessions to exclusion list...", new=new_exclusions)

        excluded_set.update(new_exclusions)

        listing["accessions"]["excluded"] = list(excluded_set)

        async with aiofiles.open(listing_path, "w") as f:
            await f.write(json.dumps(listing, indent=2, sort_keys=True))

        logger.debug("Wrote listing to file", filename=listing_path.name)


def in_refseq(comments: str, excluded: set) -> str:
    """
    Inspects text from the GenBank Comment field of a RefSeq record and retrieves
    an accession that this RefSeq record is "identical" or "derived from".

    :param comments: A string retrieved from a RefSeq record's Comment field
    :param excluded: A set of excluded listings already present in the catalog listing
    :return: If found, returns the duplicate accession, otherwise returns an empty string
    """
    for part in comments.split("."):
        if "identical" in part or "derived from" in part:
            # retrieves the last "word" from the sentence and removes the last character
            redundant_accession = part.strip().split()[-1:][0]

            if redundant_accession not in excluded:
                return redundant_accession

    return ""


async def fetch_upstream_records(
    fetch_list: list, logger: BoundLogger = base_logger
) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank

    :param fetch_list: List of accession numbers to fetch from GenBank
    :param logger: Optional entry point for an existing BoundLogger
    :return: A list of GenBank data converted from XML to dicts if possible,
        else an empty list
    """

    try:
        handle = Entrez.efetch(db="nucleotide", id=fetch_list, rettype="gb")
    except HTTPError:
        return []

    parsed = SeqIO.parse(handle, "gb")

    try:
        record_dict = SeqIO.to_dict(parsed)
    except ValueError:
        logger.error("Found two copies of the same record. Remove duplicate from this.")
        return []
    except Exception as e:
        logger.exception(e)
        return []

    if record_dict is None:
        return []

    try:
        accession_list = [record for record in record_dict.values() if record.seq]
        return accession_list
    except Exception as e:
        logger.exception(e)
        return []
