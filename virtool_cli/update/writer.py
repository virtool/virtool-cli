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
    cache_path: Path,
    dry_run: bool = False
):
    """
    Awaits new sequence data for each OTU and writes new data into JSON files with unique Virtool IDs

    :param src_path: Path to a reference directory
    :param queue: Queue holding formatted sequence and isolate data processed by this loop
    """
    logger = structlog.get_logger()
    logger.debug("Starting writer...")

    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    while True:
        packet = await queue.get()

        otu_id = packet["otu_id"]
        logger = logger.bind(otu_id=otu_id)

        sequence_data = packet["data"]

        logger.info(f"Packet {otu_id} taken from queue")

        otu_path = search_otu_by_id(otu_id, src_path)
        if not otu_path:
            logger.error("OTU path by id not found")

            await asyncio.sleep(DEFAULT_INTERVAL)
            queue.task_done()

            continue

        logger = logger.bind(otu_path=str(otu_path))

        logger.debug("Writing packet...")

        if dry_run:
            await write_summarized_update(sequence_data, otu_id, cache_path)
            cached_update_path = (cache_path / f"{otu_id}.json")
            if cached_update_path.exists():
                logger.debug(
                    "Wrote summary to cache.",
                    cached_update_path=str((cache_path / f"{otu_id}.json"))
                )
        else:
            await write_records(otu_path, sequence_data, unique_iso, unique_seq, logger)

        await asyncio.sleep(DEFAULT_INTERVAL)
        queue.task_done()


async def write_summarized_update(
    processed_updates: list, otu_id: str, cache_path: Path
):
    summary_path = cache_path / (otu_id + ".json")

    with open(summary_path, "w") as f:
        json.dump(processed_updates, f, indent=2, sort_keys=True)
