from pathlib import Path

from virtool_cli.utils.reference import get_otu_paths, get_unique_ids
from virtool_cli.utils.storage import write_records
from structlog import BoundLogger, get_logger


async def write_sequences_to_src(
    sequences: list,
    otu_path: Path,
    src_path: Path,
    logger: BoundLogger = get_logger()
):
    """

    """
    isolate_uids, sequence_uids = await get_unique_ids(get_otu_paths(src_path))

    try:
        new_sequence_paths = await write_records(
            otu_path,
            new_sequences=sequences,
            unique_iso=isolate_uids,
            unique_seq=sequence_uids,
            logger=logger,
        )

    except Exception as e:
        logger.exception(e)
        new_sequence_paths = []

    return new_sequence_paths
