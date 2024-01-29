import asyncio
from pathlib import Path

import structlog

from virtool_cli.update.update import (
    get_no_fetch_set,
    process_records,
    request_new_records,
)
from virtool_cli.update.writer import cache_new_sequences
from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.reference import get_otu_paths, get_unique_ids
from virtool_cli.utils.storage import read_otu, write_records

base_logger = structlog.get_logger()


def run(
    otu_path: Path,
    auto_evaluate: bool = False,
    dry_run: bool = False,
    debugging: bool = False,
):
    """CLI entry point for update.update.run()

    Requests updates for a single OTU directory.
    Inspects OTU contents and exclusions and requests relevant accession if found.

    :param otu_path: Path to an OTU directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param dry_run: Caching flag, writes all update data to a single file under
        "{repo_path}/.cache/updates/{otu_id}.json",
        instead of the reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)
    logger = base_logger.bind(otu_path=str(otu_path))

    if auto_evaluate:
        logger.warning(
            "Auto-evaluation is in active development and may produce false negatives.",
        )

    asyncio.run(update_otu(otu_path, auto_evaluate, dry_run))


async def update_otu(
    otu_path: Path, auto_evaluate: bool = False, dry_run: bool = False,
):
    """Requests new records for a single taxon ID
    and writes new data under the corresponding path.

    :param otu_path: Path to an OTU directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param dry_run:
    """
    src_path = otu_path.parent

    # List all isolate and sequence IDs presently in src
    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    metadata = await read_otu(otu_path)

    no_fetch_set = await get_no_fetch_set(otu_path)

    otu_id = metadata.get("_id")
    taxid = metadata.get("taxid")

    logger = base_logger.bind(taxid=taxid, otu_id=otu_id)

    try:
        record_data = await request_new_records(taxid, no_fetch_set, logger)
    except Exception as e:
        logger.exception(e)
        return

    if not record_data:
        logger.debug("No records found.")
        return

    otu_updates = await process_records(
        records=record_data,
        metadata=metadata,
        no_fetch_set=no_fetch_set,
        auto_evaluate=auto_evaluate,
        logger=logger,
    )

    if not otu_updates:
        return

    if dry_run:
        cache_path = src_path.parent / ".cache"
        update_cache_path = cache_path / "updates"
        update_cache_path.mkdir(exist_ok=True)

        await cache_new_sequences(otu_updates, otu_id, update_cache_path)

        if not update_cache_path / f"{otu_id}.json".exists():
            logger.error("Write failed")
    else:
        await write_records(otu_path, otu_updates, unique_iso, unique_seq, logger=logger)
