from pathlib import Path
import asyncio
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import (
    get_otu_paths, read_otu, 
    get_otu_accessions, get_otu_accessions_metadata
)
from virtool_cli.accessions.listings import (
    generate_listing, write_listing
)

def run(src: Path, catalog: Path, debugging: bool = False):
    """
    CLI entry point for accession.initialize.run()

    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    asyncio.run(initialize(src, catalog))

async def initialize(src_path: Path, catalog_path: Path):
    """
    Initialize a new catalog at catalog_path using data from reference directory src

    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    if not catalog_path.exists():
        catalog_path.mkdir()

    logger = base_logger.bind(src=str(src_path), catalog=str(catalog_path))

    logger.debug("Starting catalog generation...")

    queue = asyncio.Queue()

    fetcher = asyncio.create_task(
        fetcher_loop(src_path, queue))
    
    asyncio.create_task(
        writer_loop(catalog_path, queue)
    )
    
    await asyncio.gather(fetcher, return_exceptions=True)

    await queue.join()

    logger.info("Catalog generated")


async def fetcher_loop(src: Path, queue: asyncio.Queue):
    """
    Iterates through all OTUs in the src directory and generates listings for each,
    then pushes the new listing to the write queue

    :param src_path: Path to a reference directory
    :param queue: Queue holding relevant OTU information from src and fetched NCBI taxonomy id
    """
    logger = base_logger.bind(src=str(src))
    logger.debug("Starting fetcher...")
    
    for otu_path in get_otu_paths(src):
        logger = base_logger.bind(
            otu_path=str(otu_path.name)
        )

        otu_data = read_otu(otu_path)
        otu_id = otu_data['_id']
        
        logger = logger.bind(
            otu_name=otu_data.get('name', ''),
            otu_id=otu_id,
        )
        
        sequences = get_otu_accessions_metadata(otu_path)
        accessions = list(sequences.keys())
        
        listing = await generate_listing(
            otu_data=otu_data, 
            accession_list=accessions, 
            sequence_metadata=sequences, 
            logger=logger
        )

        await queue.put({ 'otu_id': otu_id, 'listing': listing } )

async def writer_loop(catalog: Path, queue: asyncio.Queue) -> None:
    """
    Pulls packet dicts from the queue and calls the write function

    :param src_path: Path to a given reference directory
    :param queue: Queue of parsed OTU data awaiting processing
    """
    while True:
        packet = await queue.get()

        logger = base_logger.bind(
            catalog=str(catalog),
            otu_id=packet['otu_id']
        )

        logger.debug(f"Got listing data for {packet['otu_id']} from the queue")

        listing = packet['listing']
        taxid = listing['taxid']

        await write_listing(
            taxid=taxid, 
            listing=listing, 
            catalog_path=catalog,
            logger=logger)

        await asyncio.sleep(0.1)
        queue.task_done()