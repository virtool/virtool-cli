from pathlib import Path
import asyncio
import structlog
from structlog import get_logger

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.utils.storage import read_otu
from virtool_cli.catalog.listings import generate_listing, write_new_listing
from virtool_cli.catalog.helpers import get_otu_accessions_metadata


def run(src_path: Path, catalog_path: Path, debugging: bool = False):
    """
    CLI entry point for catalog.initialize.run()

    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = get_logger().bind(src=str(src_path), catalog=str(catalog_path))
    logger.info("Creating new catalog in catalog path")

    asyncio.run(initialize(src_path, catalog_path))


async def initialize(src_path: Path, catalog_path: Path):
    """
    Initialize a new catalog at catalog_path using data from reference directory src

    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    if not catalog_path.exists():
        catalog_path.mkdir()

    logger = get_logger().bind(src=str(src_path), catalog=str(catalog_path))

    logger.debug("Starting catalog generation...")

    queue = asyncio.Queue()

    fetcher = asyncio.create_task(fetcher_loop(src_path, queue))

    asyncio.create_task(writer_loop(catalog_path, queue))

    await asyncio.gather(fetcher, return_exceptions=True)

    await queue.join()

    logger.info("Catalog generated")


async def fetcher_loop(src_path: Path, queue: asyncio.Queue):
    """
    Iterates through all OTUs in the src directory and generates listings for each,
    then pushes the new listing to the write queue

    :param src_path: Path to a reference directory
    :param queue: Queue holding relevant OTU information from src and fetched NCBI taxonomy id
    """
    logger = get_logger(__name__ + ".fetcher").bind(src=str(src_path))
    logger.debug("Starting fetcher...")

    for otu_path in get_otu_paths(src_path):
        logger = logger.bind(otu_path=str(otu_path.name))

        otu_data = await read_otu(otu_path)
        otu_id = otu_data["_id"]

        logger = logger.bind(
            otu_name=otu_data.get("name", ""),
            otu_id=otu_id,
        )

        sequences = get_otu_accessions_metadata(otu_path)
        accessions = list(sequences.keys())

        listing = await generate_listing(
            otu_data=otu_data,
            accession_list=accessions,
            sequence_metadata=sequences,
            logger=logger,
        )

        await queue.put({"otu_id": otu_id, "listing": listing})


async def writer_loop(catalog_path: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the write function

    :param catalog_path: Path to an accession catalog directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    logger = get_logger(__name__ + ".writer").bind(catalog=str(catalog_path))

    while True:
        packet = await queue.get()

        logger = logger.bind(otu_id=packet["otu_id"])

        logger.debug(f"Got listing data for {packet['otu_id']} from the queue")

        listing = packet["listing"]
        listing["accessions"]["excluded"] = {}

        listing_path = await write_new_listing(listing, catalog_path)

        if listing_path is None:
            logger.error("Listing could not be created under catalog")

        await asyncio.sleep(0.1)
        queue.task_done()
