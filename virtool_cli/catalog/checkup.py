from pathlib import Path
import json
import asyncio
import structlog
from structlog import get_logger
from structlog import BoundLogger

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.ncbi import get_spelling
from virtool_cli.catalog.catalog import get_catalog_paths

LISTING_KEYS = {"_id", "accessions", "name", "schema", "taxid"}


def run(catalog_path: Path, debugging: bool = False):
    """
    CLI entry point for catalog.checkup.run()

    :param catalog_path: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)

    run_tests(catalog_path)


def run_tests(catalog_path: Path):
    """
    Check catalog for outstanding issues.
    Can also be used to diagnose issues with the corresponding reference.

    :param catalog_path: Path to an accession catalog directory
    """
    logger = get_logger().bind(catalog=str(catalog_path))

    logger.info("Checking for missing data within all listings...")
    check_missing_data(catalog_path)

    logger.info("Checking for listings without taxon IDs...", test="no_taxid")
    if unassigned := search_by_taxid("none", catalog_path):
        logger.warning(
            "Found listings without assigned NCBI taxon IDs.",
            test="no_taxid",
            unassigned=unassigned,
        )

    logger.info("Checking for duplicate OTUs...", test="duplicate_taxid")
    if duplicate_otus := find_shared_taxids(catalog_path, logger):
        logger.warning(
            "Found non-unique taxon IDs in catalog.",
            test="duplicate_taxid",
            taxids=duplicate_otus,
        )

    logger.info(
        "Checking for listings containing duplicate accessions...",
        test="duplicate_accessions",
    )
    if duplicate_accessions := find_duplicate_accessions(catalog_path):
        logger.warning(
            "Found non-unique accessions in listings",
            test="duplicate_accessions",
            otus=duplicate_accessions,
        )

    logger.info(
        "Requesting spelling suggestions for unmatchable listings...",
        test="suggest_spellings",
    )
    asyncio.run(suggest_spellings(catalog_path))


def check_missing_data(catalog_path: Path):
    """
    Checks all catalog listings for missing keys and schema data.

    :param catalog_path: Path to a catalog directory
    """
    logger = get_logger(__name__ + ".missing_data")
    for listing_path in get_catalog_paths(catalog_path):
        logger = logger.bind(listing=listing_path.stem)

        listing = json.loads(listing_path.read_text())

        if missing_keys := check_keys(listing):
            logger.warning(
                "Entry is missing keys",
                missing_keys=missing_keys,
                path=str(listing_path.relative_to(catalog_path, logger)),
            )

        if not check_schema(listing):
            logger.warning("Schema field is empty.")


def check_keys(listing: dict):
    """
    Checks a listing to ensure that all necessary keys are present.
    Returns an empty list if all keys are present.

    :param listing: Catalog listing data in dictionary form
    :returns: A list of all missing keys
    """
    if not LISTING_KEYS.issubset(listing.keys()):
        missing_keys = LISTING_KEYS.difference(listing.keys())
        return list(missing_keys)

    return []


def check_schema(listing: dict):
    """
    Checks a listing to ensure that the schema field has data in it and returns a boolean

    :param listing: Catalog listing data in dictionary form
    :returns: A boolean confirming the existence of a filled-out schema
    """
    if "schema" not in listing:
        return False

    schema = listing["schema"]

    if not schema:
        return False

    if type(schema) != list:
        return False

    return True


def find_shared_taxids(catalog_path: Path, logger: BoundLogger = get_logger()) -> list:
    """
    Go through a catalog path and find listings that share a taxon ID
    and return those taxon IDs as a list

    :param catalog_path: Path to a catalog directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: A list of taxon IDs that have multiple OTUs under them
    """
    duplicated_taxids = set()

    for listing_path in catalog_path.glob("*--*.json"):
        [taxid, otu_id] = listing_path.stem.split("--")

        logger = logger.bind(
            path=str(listing_path.relative_to(catalog_path.parent)), otu_id=otu_id
        )

        if taxid != "none":
            matches = search_by_taxid(taxid, catalog_path)

            if len(matches) > 1:
                logger.debug("Duplicate taxon id found", taxid=taxid)
                duplicated_taxids.add(taxid)

    return list(duplicated_taxids)


def find_duplicate_accessions(catalog_path: Path):
    """
    Checks catalog listings for accessions that have been listed twice.
    Can be used to identify redundant sequences in the reference directory.

    :param catalog_path: Path to a catalog directory
    """
    logger = get_logger()

    duplicate_accessions = []

    for listing_path in catalog_path.glob("*--*.json"):
        listing = json.loads(listing_path.read_text())

        logger = logger.bind(listing=listing_path.name, name=listing["name"])

        logger.debug(
            "Checking for non-unique accessions within included/excluded lists..."
        )
        for alist_type in ["included", "excluded"]:
            accession_list = listing["accessions"][alist_type]
            accession_set = set(accession_list)

            if len(accession_set) < len(accession_list):
                logger.warning("Contains non-unique accessions")
                duplicate_accessions.append(listing_path.name)

        logger.debug("Checking included list against excluded list for eliminations...")
        for versioned_accession in listing["accessions"]["included"]:
            accession = versioned_accession.split(".")[0]

            if accession in accession_set:
                logger.warning(
                    f"Included accession '{versioned_accession}' is on the exclusion list."
                )
                duplicate_accessions.append(listing_path.name)

    return duplicate_accessions


async def suggest_spellings(catalog_path: Path):
    """
    Evaluates the names of OTUs without retrievable taxon IDs
    and queries Entrez ESpell for alternatives.

    :param catalog_path: Path to a catalog directory
    """
    logger = get_logger()

    for listing_path in catalog_path.glob("none--*.json"):
        with open(listing_path, "r") as f:
            listing = json.load(f)

        logger = logger.bind(
            path=str(listing_path.relative_to(catalog_path)),
            otu_id=listing["_id"],
            taxid=listing["taxid"],
            current_name=listing["name"],
        )

        current_name = listing["name"]

        if "-" in current_name:
            current_name = current_name.split("-")[0]

        new_spelling = await get_spelling(current_name)

        if new_spelling != current_name.lower():
            logger.info(f"Try: {new_spelling}", potential_name=new_spelling)


def search_by_taxid(taxid: int | str, catalog_path: Path) -> list:
    """
    Searches records for a matching taxon id and returns all matching paths in the accession records as strings
    (for logging purposes)

    :param taxid: Unique taxon ID
    :param catalog_path: Path to an accession catalog directory
    :return: List of listings matching this taxon ID
    """
    matches = [
        str(listing.relative_to(catalog_path))
        for listing in catalog_path.glob(f"{taxid}--*.json")
    ]

    if matches:
        return matches

    return []
