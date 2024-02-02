from pathlib import Path

import structlog

from virtool_cli.check.otu import check_otu
from virtool_cli.utils.reference import get_otu_paths

base_logger = structlog.get_logger()


def check_reference(
    src_path: Path,
    logger: structlog.BoundLogger = base_logger,
) -> bool:
    """Inspects all files in the reference directory and logs a warning if a problem is found.

    :param src_path: Path to a given reference directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    ref_valid = True

    logger.info("Running checks on reference directory...")

    watch_list = []

    for otu_path in get_otu_paths(src_path):
        otu_valid = check_otu(otu_path)

        if not otu_valid:
            ref_valid = False
            watch_list.append(otu_path.name)

    # Print outstanding issues to stdout
    print(watch_list)

    return ref_valid
