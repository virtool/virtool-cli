from pathlib import Path
import structlog
from virtool_cli.utils.logging import DEBUG_LOGGER, DEFAULT_LOGGER
from virtool_cli.check.checkup import check_otu

base_logger = structlog.get_logger()


def run(otu_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.checkup.run()

    :param otu_path: Path to an OTU directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(otu_path=str(otu_path), verbose=debugging)

    logger.debug("Debug flag is enabled")

    if check_otu(otu_path):
        logger.info("Reference folder is valid.", valid_result=True)
        print(True)
    else:
        logger.error("Reference folder is invalid.", valid_result=False)
        print(False)
