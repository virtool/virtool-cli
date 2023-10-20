from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import get_otu_paths, get_unique_ids
from virtool_cli.utils.format import process_records
from virtool_cli.utils.storage import write_records
from virtool_cli.catalog.listings import parse_listing
from virtool_cli.catalog.helpers import search_by_otu_id
from virtool_cli.update.update import request_new_records

base_logger = structlog.get_logger()


def run(
    otu_path: Path,
    catalog_path: Path,
    auto_evaluate: bool = False,
    debugging: bool = False,
):
    """
    CLI entry point for update.update.run()

    Requests updates for a single OTU directory
    Searches the catalog for a matching catalog listing and requests new updates if it finds one.

    :param otu_path: Path to a OTU directory
    :param catalog_path: Path to a catalog directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(otu_path=str(otu_path))

    if auto_evaluate:
        logger.warning(
            "Auto-evaluation is in active development and may produce false negatives."
        )

    [_, otu_id] = otu_path.name.split("--")
    listing_path = search_by_otu_id(otu_id, catalog_path)

    if listing_path:
        asyncio.run(update_otu(otu_path, listing_path, auto_evaluate))

    else:
        logger.error("Listing not found for this OTU.")
        logger.info("Use virtool acc init to create a new accession catalog.")


async def update_otu(otu_path: Path, listing_path: Path, auto_evaluate: bool = False):
    """
    Requests new records for a single taxon ID
    and writes new data under the corresponding path.

    :param otu_path: Path to a OTU directory
    :param listing_path: Path to a listing in an accession catalog directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    """
    src_path = otu_path.parent
    listing = await parse_listing(listing_path)

    # extract taxon ID and _id hash from listing filename
    [taxid, otu_id] = listing_path.stem.split("--")

    logger = base_logger.bind(taxid=taxid, otu_id=otu_id)

    record_data = await request_new_records(listing)
    if not record_data:
        return

    otu_updates = await process_records(
        records=record_data, listing=listing, auto_evaluate=auto_evaluate, logger=logger
    )
    if not otu_updates:
        return

    # List all isolate and sequence IDs presently in src
    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    await write_records(otu_path, otu_updates, unique_iso, unique_seq, logger=logger)
