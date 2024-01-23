import json
from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.reference import (
    get_otu_paths,
    search_otu_by_id,
    get_unique_ids,
)
from virtool_cli.utils.storage import write_records

DEFAULT_INTERVAL = 0.001


async def writer_loop(
    src_path: Path,
    queue: asyncio.Queue,
):
    """
    Awaits new sequence data per OTU and writes new data
    to the correct location under the reference directory

    :param src_path: Path to a reference directory
    :param queue: Queue holding formatted sequence and isolate data processed by this loop
    """
    logger = structlog.get_logger()
    logger.debug("Starting writer...")

    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    while True:
        packet = await queue.get()
        otu_id, sequence_data = await process_packet(packet)

        logger = logger.bind(otu_id=otu_id)

        sequence_data = packet["data"]

        otu_path = await get_otu_path(otu_id, src_path, logger)
        if not otu_path:
            queue.task_done()
            continue

        logger = logger.bind(otu_path=str(otu_path))

        logger.debug("Writing packet...")
        await write_records(
            otu_path=otu_path,
            new_sequences=sequence_data,
            unique_iso=unique_iso,
            unique_seq=unique_seq,
            logger=logger
        )

        await asyncio.sleep(DEFAULT_INTERVAL)
        queue.task_done()


async def cacher_loop(
    src_path: Path,
    cache_path: Path,
    queue: asyncio.Queue,
):
    """
    Awaits new sequence data per OTU and writes new data into one JSON file per OTU

    :param src_path: Path to a reference directory
    :param cache_path: Path to a directory containing cached update lists
    :param queue: Queue holding formatted sequence and isolate data processed by this loop
    """
    logger = structlog.get_logger()
    logger.debug("Starting cache writer...")

    while True:
        packet = await queue.get()
        otu_id, sequence_data = await process_packet(packet)

        logger = logger.bind(otu_id=otu_id)

        logger.info(f"Packet {otu_id} taken from queue")

        otu_path = await get_otu_path(otu_id, src_path, logger)
        if not otu_path:
            queue.task_done()
            continue

        logger = logger.bind(otu_path=str(otu_path))
        logger.debug("Writing packet...")

        await cache_new_sequences(sequence_data, otu_id, cache_path)
        cached_update_path = (cache_path / f"{otu_id}.json")
        if cached_update_path.exists():
            logger.debug(
                "Wrote summary to cache.",
                cached_update_path=str((cache_path / f"{otu_id}.json"))
            )

        await asyncio.sleep(DEFAULT_INTERVAL)
        queue.task_done()


async def get_otu_path(otu_id, src_path, logger):
    otu_path = search_otu_by_id(otu_id, src_path)
    if not otu_path:
        logger.error("OTU path by id not found")

        await asyncio.sleep(DEFAULT_INTERVAL)
        return None

    return otu_path

async def process_packet(packet):
    otu_id = packet["otu_id"]
    sequence_data = packet["data"]

    return otu_id, sequence_data


async def cache_new_sequences(
    processed_updates: list[dict], otu_id: str, cache_path: Path
):
    """
    Takes a list of processed update data and caches it under a given cache path

    :param processed_updates: Preprocessed sequence records
    :param otu_id: Unique OTU identifier
    :param cache_path: Path to a directory containing cached update lists
    """
    summary_path = cache_path / (otu_id + ".json")

    with open(summary_path, "w") as f:
        json.dump(processed_updates, f, indent=2, sort_keys=True)
