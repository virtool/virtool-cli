from pathlib import Path
import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.doctor.fix_otu import repair_otu

base_logger = structlog.get_logger()


def run(src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.fix_reference.run()

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(src=str(src_path))

    logger.info("Repairing data...")

    repair_data(src_path)


def repair_data(src_path: Path, logger: structlog.BoundLogger = base_logger):
    """
    Fixes incorrect data in all OTUs under a reference

    :param src_path: Path to a given reference directory
    :param logger: Optional entry point for an existing BoundLogger
    """
    for otu_path in get_otu_paths(src_path):
        repair_otu(otu_path, logger)
