from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import (
    get_otu_paths,
    search_otu_by_id,
    is_v1,
    get_unique_ids,
)
from virtool_cli.utils.format import process_records
from virtool_cli.utils.storage import write_records, read_otu
from virtool_cli.update.update import request_new_records, get_no_fetch_set

DEFAULT_INTERVAL = 0.001

base_logger = structlog.get_logger()


def run(
    src_path: Path,
    auto_evaluate: bool = False,
    debugging: bool = False,
):
    """
    CLI entry point for update.update_ref.run()

    Requests updates for all OTU directories under a source reference

    :param src_path: Path to a reference directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
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
        update_reference(src_path=src_path, auto_evaluate=auto_evaluate)
    )


async def update_reference(
    src_path: Path, auto_evaluate: bool = False
):
    """
    Creates 2 queues:
        1) upstream: Holds raw NCBI GenBank data,
        2) write: Holds formatted sequence data and isolate data

    Monitors 3 asynchronous processes:
        1) fetcher:
            Requests and retrieves new accessions from NCBI GenBank
            and pushes results to upstream queue
        2) processor:
            Pulls Genbank data from upstream queue,
            formats into dict form and pushes to write queue
        3) writer:
            Pulls formatted sequences from write queue
            and writes json to the correct location in the
            src directory

    :param src_path: Path to a reference directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    """
    # Holds raw NCBI GenBank data
    upstream_queue = asyncio.Queue()

    # Holds formatted sequence data and isolate data
    write_queue = asyncio.Queue()

    # Requests and retrieves new accessions from NCBI GenBank
    # and pushes results to upstream queue
    fetcher = asyncio.create_task(
        fetcher_loop(get_otu_paths(src_path), queue=upstream_queue)
    )

    # Pulls Genbank data from upstream queue, formats into dict form
    # and pushes to write queue
    asyncio.create_task(
        processor_loop(upstream_queue, write_queue, auto_evaluate=auto_evaluate)
    )

    # Pulls formatted sequences from write queue, checks isolate metadata
    # and writes json to the correct location in the src directory
    asyncio.create_task(writer_loop(src_path, write_queue))

    await asyncio.gather(*[fetcher], return_exceptions=True)

    await upstream_queue.join()  # wait until the consumer has processed all items

    await write_queue.join()

    return


async def fetcher_loop(otu_paths: list, queue: asyncio.Queue):
    """
    Loops through selected OTU listings from accession catalogue,
    indexed by NCBI taxon ID, and:
        1) requests NCBI Genbank for accession numbers not extant
            in catalog,
        2) loops through retrieved new accession numbers and
            requests relevant record data from NCBI Genbank
        3) Pushes new records and corresponding OTU information
            to a queue for formatting

    :param otu_paths: A list of OTU paths
    :param queue: Queue holding fetched NCBI GenBank data
    """
    logger = structlog.get_logger(__name__ + ".fetcher")
    logger.debug("Starting fetcher...")

    for path in otu_paths:
        otu_metadata = await read_otu(path)

        otu_id = otu_metadata.get('_id')
        taxid = otu_metadata.get('taxid')

        logger = logger.bind(taxid=taxid, otu_id=otu_id)
        logger.debug("Starting OTU...")

        try:
            no_fetch_set = await get_no_fetch_set(path)
        except Exception as e:
            logger.exception(e)
            continue

        record_data = await request_new_records(taxid, no_fetch_set, logger)
        if not record_data:
            continue

        packet = {
            "taxid": taxid,
            "otu_id": otu_id,
            "metadata": otu_metadata,
            "no_fetch_set": no_fetch_set,
            "listing": otu_metadata,
            "data": record_data,
        }

        await queue.put(packet)
        logger.debug(
            f"Pushed {len(record_data)} requests to upstream queue",
            n_requests=len(record_data),
            taxid=taxid,
        )
        await asyncio.sleep(DEFAULT_INTERVAL)


async def processor_loop(
    upstream_queue: asyncio.Queue,
    downstream_queue: asyncio.Queue,
    auto_evaluate: bool = False,
):
    """
    Awaits fetched sequence data from the fetcher:
        1) Formats the sequence data into reference-compatible dictionaries,
        2) Checks for validity,
        3) Pushes the formatted data into the downstream queue to be dealt with by the writer

    :param upstream_queue: Queue holding NCBI GenBank data pushed by the fetcher,
    :param downstream_queue: Queue holding formatted sequence and isolate data processed by this loop
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    """
    logger = structlog.get_logger(__name__ + ".processor")
    logger.debug("Starting processor...")

    while True:
        fetch_packet = await upstream_queue.get()

        metadata = fetch_packet["metadata"]
        otu_id = metadata["_id"]
        taxid = metadata.get("taxid")
        logger = logger.bind(taxid=taxid)

        otu_updates = await process_records(
            records=fetch_packet["data"],
            metadata=metadata,
            no_fetch_set=fetch_packet["no_fetch_set"],
            auto_evaluate=auto_evaluate,
            logger=logger,
        )

        if not otu_updates:
            await asyncio.sleep(DEFAULT_INTERVAL)
            upstream_queue.task_done()
            continue

        processed_packet = {"taxid": taxid, "otu_id": otu_id, "data": otu_updates}

        await downstream_queue.put(processed_packet)
        logger.debug(f"Pushed {len(otu_updates)} new accessions to downstream queue")
        await asyncio.sleep(DEFAULT_INTERVAL)
        upstream_queue.task_done()


async def writer_loop(
    src_path: Path,
    queue: asyncio.Queue,
):
    """
    Awaits new sequence data for each OTU and writes new data into JSON files with unique Virtool IDs

    :param src_path: Path to a reference directory
    :param queue: Queue holding formatted sequence and isolate data processed by this loop
    """
    logger = structlog.get_logger(__name__ + ".writer")
    logger.debug("Starting writer...")

    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    while True:
        packet = await queue.get()

        taxid = packet["taxid"]
        otu_id = packet["otu_id"]
        sequence_data = packet["data"]

        logger = logger.bind(otu_id=otu_id, taxid=taxid, src=str(src_path))

        otu_path = search_otu_by_id(otu_id, src_path)

        await write_records(otu_path, sequence_data, unique_iso, unique_seq, logger)

        await asyncio.sleep(DEFAULT_INTERVAL)
        queue.task_done()
