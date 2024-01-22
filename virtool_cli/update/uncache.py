from pathlib import Path
import json
import asyncio
import structlog

from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.reference import is_v1
from virtool_cli.update.writer import writer_loop

DEFAULT_INTERVAL = 0.001

base_logger = structlog.get_logger()


def run(
    cache_path: Path,
    src_path: Path,
    auto_evaluate: bool = False,
    debugging: bool = False,
):
    """
    CLI entry point for update.update_ref.run()

    Requests updates for all OTU directories under a source reference

    :param src_path: Path to a reference directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param dry_run:
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)
    logger = base_logger.bind(src=str(src_path))

    if is_v1(src_path):
        logger.error(
            'reference folder "src" is a deprecated v1 reference.'
            + 'Run "virtool ref migrate" before trying again.'
        )
        return

    if auto_evaluate:
        logger.warning(
            "Auto-evaluation is in active development and may produce false negatives."
        )

    logger.info("Updating src directory accessions...")

    asyncio.run(
        update_reference_from_cache(
            cache_path=cache_path, src_path=src_path, auto_evaluate=auto_evaluate
        )
    )


async def update_reference_from_cache(
    cache_path: Path, src_path: Path, auto_evaluate: bool = False
):
    """
    """
    # Holds raw NCBI GenBank data
    queue = asyncio.Queue()

    # Requests and retrieves new accessions from NCBI GenBank
    # and pushes results to upstream queue
    fetcher = asyncio.create_task(
        loader_loop(cache_path=cache_path, queue=queue)
    )

    # Pulls formatted sequences from write queue, checks isolate metadata
    # and writes json to the correct location in the src directory
    asyncio.create_task(writer_loop(src_path, queue, cache_path, dry_run=False))

    await asyncio.gather(*[fetcher], return_exceptions=True)

    await queue.join()

    return


async def loader_loop(
    cache_path: Path, queue: asyncio.Queue
):
    """
    Loops through selected OTU listings from accession catalogue,
    indexed by NCBI taxon ID, and:
        1) requests NCBI Genbank for accession numbers not extant
            in catalog,
        2) loops through retrieved new accession numbers and
            requests relevant record data from NCBI Genbank
        3) Pushes new records and corresponding OTU information
            to a queue for formatting

    :param cache_path:
    :param queue: Queue holding fetched NCBI GenBank data
    """
    logger = structlog.get_logger(__name__ + ".fetcher")
    logger.debug("Starting loader...")

    for path in cache_path.glob("*.json"):
        otu_id = path.stem

        logger = logger.bind(otu_id = otu_id)
        logger.debug(f"Loading {otu_id}...")

        with open(path, "r") as f:
            record_data = json.load(f)

        packet = {
            "otu_id": otu_id,
            "data": record_data,
        }

        await queue.put(packet)
        logger.debug(
            f"Pushed {len(record_data)} requests to queue",
            n_requests=len(record_data)
        )
        await asyncio.sleep(DEFAULT_INTERVAL)