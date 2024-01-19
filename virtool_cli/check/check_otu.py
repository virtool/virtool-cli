from pathlib import Path
import structlog

from virtool_cli.utils.logging import configure_logger
from virtool_cli.check.checkup import check_otu

base_logger = structlog.get_logger()


def run(otu_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.checkup.run()

    :param otu_path: Path to an OTU directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)
    logger = base_logger.bind(otu_path=str(otu_path), verbose=debugging)

    logger.debug("Debug flag is enabled")

    if check_otu(otu_path):
        logger.info("Reference folder is valid.", valid_result=True)
        print(True)
    else:
        logger.error("Reference folder is invalid.", valid_result=False)
        print(False)
