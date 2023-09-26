from pathlib import Path
import json
import asyncio
import structlog
from structlog.stdlib import get_logger, BoundLogger

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import (
    get_otu_paths,
    search_otu_by_id,
    read_otu,
)
from virtool_cli.catalog.listings import generate_listing, write_listing
from virtool_cli.catalog.helpers import (
    get_catalog_paths,
    filter_catalog,
    get_otu_accessions,
    get_otu_accessions_metadata,
)


def run(src_path: Path, catalog_path: Path, debugging: bool = False):
    """
    CLI entry point for catalog.update.run()

    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = get_logger().bind(catalog=str(catalog_path))

    logger.info("Starting catalog updater...")

    asyncio.run(update(src_path, catalog_path))


async def update(src_path: Path, catalog_path: Path):
    """
    Checks reference directory src for data not accounted for by the accession catalog
    and updates the corresponding catalog listing accordingly.

    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    logger = get_logger().bind(task="update", catalog=str(catalog_path))

    await complete_catalog(src_path, catalog_path, logger)

    queue = asyncio.Queue()

    fetcher = asyncio.create_task(fetcher_loop(src_path, catalog_path, queue))

    asyncio.create_task(writer_loop(catalog_path, queue))

    await asyncio.gather(fetcher, return_exceptions=True)

    await queue.join()

    logger.info("Catalog up to date")


async def fetcher_loop(src_path: Path, catalog_path: Path, queue: asyncio.Queue):
    """
    Iterates through all listings with an associated OTU directory and:
        1) Inspects the OTU directory for changes.
        2) Writes changes to the catalog listing
        3) Pushes the changes to the queue for writing

    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    :param queue: Queue containing changed listings
    """
    fetch_logger = get_logger().bind(loop="fetcher")
    fetch_logger.debug("Starting fetcher...")

    update_count = 0

    for listing_path in filter_catalog(src_path, catalog_path):
        [taxid, otu_id] = listing_path.stem.split("--")

        logger = fetch_logger.bind(listing=listing_path.name)

        otu_path = search_otu_by_id(otu_id, src_path)
        existing_accessions = set(get_otu_accessions(otu_path))
        fetch_logger.debug(f"current accessions: {existing_accessions}")

        listing = json.loads(listing_path.read_text())

        logger = logger.bind(otu_name=listing["name"], otu_id=otu_id, taxid=taxid)

        listed_accessions = set(listing["accessions"]["included"])

        not_in_reference = listed_accessions.difference(existing_accessions)
        not_in_listing = existing_accessions.difference(listed_accessions)

        indexed_accessions = listing["accessions"]["included"]

        if not_in_reference:
            logger.debug(f"Reference is missing accessions: {not_in_reference}")
            for accession in not_in_reference:
                indexed_accessions.pop(accession)

        if not_in_listing:
            logger.debug(f"Listing is missing accessions: {not_in_listing}")

            for accession in not_in_listing:
                if accession in indexed_accessions:
                    logger.debug(f"Removing missing accession {accession}")
                    indexed_accessions.pop(accession)
                else:
                    indexed_accessions.append(accession)

                listing["accessions"]["included"] = indexed_accessions

            await queue.put({"path": listing_path, "listing": listing})

            logger.debug(
                "Pushed updated listing to queue", updated_accessions=indexed_accessions
            )

            update_count += 1

    fetch_logger.debug(
        f"Pushed {update_count} updates to writer",
        n_updated=update_count,
        catalog_path=str(catalog_path),
    )


async def writer_loop(catalog_path: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the update function

    :param catalog_path: Path to a catalog directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    write_logger = get_logger().bind(loop="writer", catalog=str(catalog_path))
    write_logger.debug("Starting writer...")

    while True:
        packet = await queue.get()
        listing_path = packet["path"]
        listing = packet["listing"]
        write_logger.bind(
            listing_path=str(listing_path.relative_to(catalog_path)),
            logger=write_logger,
        )

        write_logger.debug(f"New listing:\n{listing}")

        with open(packet["path"], "w") as f:
            json.dump(packet["listing"], f, indent=2, sort_keys=True)

        await asyncio.sleep(0.1)
        queue.task_done()


async def complete_catalog(
    src_path: Path, catalog_path: Path, logger: BoundLogger = get_logger()
):
    """
    Check reference directory src contents against
    accession catalog contents for missing OTUs.
    If the reference contains an OTU that is not in the catalog,
    adds a new listing to the catalog.

    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    :param logger: Optional entry point for an existing BoundLogger
    """
    try:
        otu_ids = {otu_path.name.split("--")[1] for otu_path in get_otu_paths(src_path)}

        catalog_ids = {
            listing_path.stem.split("--")[1]
            for listing_path in get_catalog_paths(catalog_path)
        }

    except Exception as e:
        logger.exception(e)
        return

    if otu_ids.issubset(catalog_ids):
        logger.info("All OTUs present in catalog")

    # If if we're missing OTUs, keep going
    missing_otus = otu_ids.difference(catalog_ids)
    logger.warning("Missing listings for OTUs", missing_otus=missing_otus)

    for otu_id in missing_otus:
        otu_path = [path for path in src_path.glob(f"*--{otu_id}") if path.is_dir()][0]
        await add_listing(otu_path, catalog_path, logger=logger)


async def add_listing(
    otu_path: Path, catalog_path: Path, logger: BoundLogger = get_logger()
):
    """
    Generate and write a new listing to the catalog based on a given OTU

    :param otu_path: Path to an OTU directory in a reference directory
    :param catalog_path: Path to a catalog directory
    :param logger: Optional entry point for an existing BoundLogger
    """
    otu_data = read_otu(otu_path)

    logger.info(f"No accession record for {otu_data['taxid']}.")

    sequences = get_otu_accessions_metadata(otu_path)
    accessions = list(sequences.keys())

    new_listing = await generate_listing(
        otu_data=otu_data,
        accession_list=accessions,
        sequence_metadata=sequences,
        logger=logger,
    )

    if not new_listing:
        logger.error("Could not generate a listing for this OTU.")
        return

    await write_listing(
        otu_data["taxid"], new_listing, catalog_path=catalog_path, logger=logger
    )

    return
