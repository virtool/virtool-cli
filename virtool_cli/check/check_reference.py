from pathlib import Path
import structlog

from virtool_cli.check.checkup import check_otu
from virtool_cli.utils.logging import DEBUG_LOGGER, DEFAULT_LOGGER
from virtool_cli.utils.reference import get_otu_paths

base_logger = structlog.get_logger()


def run(src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.checkup.run()

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(src=str(src_path), verbose=debugging)

    logger.debug("Debug flag is enabled")

    if check_reference(src_path):
        logger.info("Reference folder is valid.", valid_result=True)
    else:
        logger.error("Reference folder is invalid.", valid_result=False)


def check_reference(
    src_path: Path, logger: structlog.BoundLogger = base_logger
) -> bool:
    """
    Inspects all files in the reference directory and logs a warning if a problem is found.

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
