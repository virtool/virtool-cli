from pathlib import Path
import structlog
from structlog import get_logger

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.doctor.fix_otu import repair_otu


def run(src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.fix_reference.run()

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = get_logger().bind(src=str(src_path))

    logger.info("Repairing data...")

    repair_data(src_path)


def repair_data(src_path: Path):
    """
    Fixes incorrect data in all OTUs under a reference

    :param src_path: Path to a given reference directory
    """
    for otu_path in get_otu_paths(src_path):
        repair_otu(otu_path, get_logger().bind(otu=otu_path.relative_to(src_path)))
