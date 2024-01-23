from pathlib import Path
import asyncio
import structlog
from urllib.error import HTTPError

from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.reference import is_v1
from virtool_cli.utils.storage import read_otu
from virtool_cli.update.update import request_new_records, get_no_fetch_set, process_records
from virtool_cli.update.writer import writer_loop, cacher_loop

DEFAULT_INTERVAL = 0.001

base_logger = structlog.get_logger()


def run(
    src_path: Path,
    filter: str = "*",
    auto_evaluate: bool = False,
    dry_run: bool = False,
    debugging: bool = False,
):
    """
    CLI entry point for update.update_ref.run()

    Requests updates for all OTU directories under a source reference

    :param src_path: Path to a reference directory
    :param filter: Changing this parameter defines which paths are added to otu_paths;
        e.g. `a*` will update only OTUs that start with `a`.
        defaults to `*` (wildcard)
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param dry_run: Caching flag, writes update data for each OTU
        under "../.cache/updates/{otu_id}.json",
        instead of the reference directory
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
        update_reference(
            src_path=src_path, filter=filter, auto_evaluate=auto_evaluate, dry_run=dry_run)
    )


async def update_reference(
    src_path: Path,
    filter: str = "*",
    auto_evaluate: bool = False,
    dry_run: bool = False
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
    :param filter: Filter criteria for updates
    :param auto_evaluate: Auto-evaluation flag,
        enables automatic filtering for fetched results
    :param dry_run: Caching flag, writes update data for each OTU
        under "../.cache/updates/{otu_id}.json",
        instead of the reference directory
    """
    # Create cache if necessary
    cache_path = src_path.parent / ".cache"
    update_cache_path = cache_path / "updates"
    update_cache_path.mkdir(exist_ok=True)

    # Holds raw NCBI GenBank data
    upstream_queue = asyncio.Queue()

    # Holds formatted sequence data and isolate data
    write_queue = asyncio.Queue()

    # Applies a glob filter to OTU:
    otu_paths = await filter_otu_paths(src_path, filter)

    # Requests and retrieves new accessions from NCBI GenBank
    # and pushes results to upstream queue
    fetcher = asyncio.create_task(
        fetcher_loop(
            otu_paths, queue=upstream_queue, cache_path=update_cache_path, dry_run=dry_run
        )
    )

    # Pulls Genbank data from upstream queue, formats into dict form
    # and pushes to write queue
    asyncio.create_task(
        processor_loop(upstream_queue, write_queue, auto_evaluate=auto_evaluate)
    )

    if dry_run:
        # Pulls formatted sequence list from write queue and
        # writes the list to "{otu_id}.json" under the cache path
        asyncio.create_task(
            cacher_loop(src_path, update_cache_path, write_queue)
        )
    else:
        # Pulls formatted sequences from write queue, checks isolate metadata
        # and writes json to the correct location in the src directory
        asyncio.create_task(
            writer_loop(src_path, write_queue)
        )

    await asyncio.gather(*[fetcher], return_exceptions=True)

    await upstream_queue.join()  # wait until the consumer has processed all items

    await write_queue.join()

    return


async def fetcher_loop(
    otu_paths: list, queue: asyncio.Queue, cache_path: Path, dry_run: bool = False
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

    :param otu_paths: A list of OTU paths
    :param queue: Queue holding fetched NCBI GenBank data
    :param cache_path:
    :param dry_run:
    """
    logger = structlog.get_logger(__name__ + ".fetcher")
    logger.debug("Starting fetcher...")

    for path in otu_paths:
        otu_metadata = await read_otu(path)

        logger.debug(otu_metadata)

        otu_id = otu_metadata['_id']
        taxid = otu_metadata.get('taxid', None)
        logger = logger.bind(taxid=taxid, otu_id=otu_id)

        if taxid is None:
            logger.error("NCBI Taxonomy id not found in OTU metadata. Moving on...")
            continue

        if dry_run:
            if (cache_path / f"{otu_id}.json").exists():
                logger.warning("OTU updates have already been cached. Moving on...")
                continue

        logger.debug("Starting OTU...")

        try:
            no_fetch_set = await get_no_fetch_set(path)
        except Exception as e:
            logger.exception(e)
            continue

        try:
            record_data = await request_new_records(taxid, no_fetch_set, logger)
        except HTTPError as e:
            logger.error(f"{e}: Failed to retrieve new records.")
            logger.exception(e)
            continue

        if not record_data:
            logger.debug("No records found.")
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

        processed_packet = {"otu_id": otu_id, "data": otu_updates}

        await downstream_queue.put(processed_packet)
        logger.debug(f"Pushed {len(otu_updates)} new accessions to downstream queue")

        await asyncio.sleep(DEFAULT_INTERVAL)
        upstream_queue.task_done()


async def filter_otu_paths(src_path: Path, filter: str = "*") -> list[Path]:
    """
    Takes a glob-formatted filter on directory names

    :param src_path:
    :param filter:
    """
    filtered_paths = []
    for path in src_path.glob(f"{filter}--*"):
        if path.is_dir():
            filtered_paths.append(path)

    return filtered_paths

