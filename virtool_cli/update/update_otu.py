from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import get_otu_paths, get_unique_ids
from virtool_cli.utils.format import process_records
from virtool_cli.utils.storage import read_otu, write_records
from virtool_cli.update.update import get_no_fetch_set, request_new_records

base_logger = structlog.get_logger()


def run(
    otu_path: Path,
    auto_evaluate: bool = False,
    debugging: bool = False,
):
    """
    CLI entry point for update.update.run()

    Requests updates for a single OTU directory
    Searches the catalog for a matching catalog listing and requests new updates if it finds one.

    :param otu_path: Path to an OTU directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(otu_path=str(otu_path))

    if auto_evaluate:
        logger.warning(
            "Auto-evaluation is in active development and may produce false negatives."
        )

    asyncio.run(update_otu(otu_path, auto_evaluate))


async def update_otu(
    otu_path: Path, auto_evaluate: bool = False
):
    """
    Requests new records for a single taxon ID
    and writes new data under the corresponding path.

    :param otu_path: Path to a OTU directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    """
    src_path = otu_path.parent
    metadata = await read_otu(otu_path)

    no_fetch_set = await get_no_fetch_set(otu_path)

    otu_id = metadata.get('_id')
    taxid = metadata.get('taxid')

    logger = base_logger.bind(taxid=taxid, otu_id=otu_id)

    try:
        record_data = await request_new_records(taxid, no_fetch_set, logger)
    except Exception as e:
        logger.exception(e)
        return

    if not record_data:
        return

    otu_updates = await process_records(
        records=record_data,
        metadata=metadata,
        no_fetch_set=no_fetch_set,
        auto_evaluate=auto_evaluate,
        logger=logger
    )
    if not otu_updates:
        return

    # List all isolate and sequence IDs presently in src
    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    await write_records(otu_path, otu_updates, unique_iso, unique_seq, logger=logger)
